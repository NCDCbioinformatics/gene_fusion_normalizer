import difflib, re

PRIORITY = ['name-exact','synonym-exact','special-intergenic','special-lnc','grouped-semicolon','fuzzy-name','fuzzy-synonym','synonym-only','unmatched']
PRIORITY_RANK = {k:i for i,k in enumerate(PRIORITY)}

SPECIAL_LOC_RE = re.compile(r'^(LOC\d+|LINC\d+|MIR\d+[A-Z]?\d*|SNOR[A-Z]+\d+.*)$', re.IGNORECASE)

class GeneResolver:
    def __init__(self, gene_name_to_ids, gene_id_to_name, synonym_map=None):
        self.gene_name_to_ids = gene_name_to_ids
        self.gene_id_to_name = gene_id_to_name
        self.synonym_map = synonym_map or {}
        self._universe = set(gene_name_to_ids.keys()) | set(self.synonym_map.keys())

    def is_grouped_semicolon(self, token: str) -> bool:
        return ';' in token

    def is_intergenic(self, token: str) -> bool:
        return token.strip().lower() == 'intergenic'

    def is_special_lnc(self, token: str) -> bool:
        return bool(SPECIAL_LOC_RE.match(token.strip()))

    def resolve_one(self, token: str, fuzzy=True, cutoff=0.86):
        if token is None:
            return {'token':'', 'matched_symbol':'', 'ensembl_gene_id':'', 'match_type':'unmatched'}
        raw = token.strip()
        key = raw.lower()
        if not key:
            return {'token':token, 'matched_symbol':'', 'ensembl_gene_id':'', 'match_type':'unmatched'}

        if self.is_grouped_semicolon(raw):
            return {'token': token, 'matched_symbol': raw, 'ensembl_gene_id': '', 'match_type': 'grouped-semicolon'}

        if self.is_intergenic(raw):
            return {'token': token, 'matched_symbol': 'INTERGENIC', 'ensembl_gene_id': '', 'match_type': 'special-intergenic'}

        if self.is_special_lnc(raw):
            if key in self.gene_name_to_ids:
                ids = self.gene_name_to_ids[key]
                gid_list = sorted(ids)
                proper = self.gene_id_to_name[gid_list[0]]
                return {'token':token,'matched_symbol':proper,'ensembl_gene_id':';'.join(gid_list),'match_type':'name-exact'}
            return {'token': token, 'matched_symbol': raw, 'ensembl_gene_id': '', 'match_type': 'special-lnc'}

        if key in self.synonym_map:
            canonical = self.synonym_map[key]
            ids = self.gene_name_to_ids.get(canonical.lower(), set())
            if ids:
                gid_list = sorted(ids)
                return {'token':token,'matched_symbol':canonical,'ensembl_gene_id':';'.join(gid_list),'match_type':'synonym-exact'}
            return {'token':token,'matched_symbol':canonical,'ensembl_gene_id':'','match_type':'synonym-only'}

        if key in self.gene_name_to_ids:
            ids = self.gene_name_to_ids[key]
            gid_list = sorted(ids)
            proper = self.gene_id_to_name[gid_list[0]]
            return {'token':token,'matched_symbol':proper,'ensembl_gene_id':';'.join(gid_list),'match_type':'name-exact'}

        if fuzzy and self._universe:
            cand = difflib.get_close_matches(key, list(self._universe), n=1, cutoff=cutoff)
            if cand:
                ck = cand[0]
                if ck in self.synonym_map:
                    canonical = self.synonym_map[ck]
                    ids = self.gene_name_to_ids.get(canonical.lower(), set())
                    gid_list = sorted(ids)
                    return {'token':token,'matched_symbol':canonical,'ensembl_gene_id':';'.join(gid_list),'match_type':'fuzzy-synonym'}
                if ck in self.gene_name_to_ids:
                    ids = self.gene_name_to_ids[ck]
                    gid_list = sorted(ids)
                    proper = self.gene_id_to_name[gid_list[0]]
                    return {'token':token,'matched_symbol':proper,'ensembl_gene_id':';'.join(gid_list),'match_type':'fuzzy-name'}

        return {'token':token,'matched_symbol':'','ensembl_gene_id':'','match_type':'unmatched'}

    @staticmethod
    def best_of(results):
        if not results:
            return None
        return sorted(results, key=lambda r: PRIORITY_RANK.get(r.get('match_type','unmatched'), 999))[0]