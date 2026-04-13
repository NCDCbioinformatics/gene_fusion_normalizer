import gzip, io, re, pandas as pd

def open_maybe_gzip(path: str):
    if path.endswith('.gz'):
        return io.TextIOWrapper(gzip.open(path, 'rb'))
    return open(path, 'r', encoding='utf-8', errors='ignore')

def parse_gtf_gene_map(gtf_path: str):
    gene_name_to_ids = {}
    gene_id_to_name = {}
    with open_maybe_gzip(gtf_path) as fh:
        for line in fh:
            if not line or line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9 or parts[2] != 'gene':
                continue
            attrs = parts[8]
            m_id = re.search(r'gene_id "([^"]+)"', attrs)
            m_name = re.search(r'gene_name "([^"]+)"', attrs)
            if not m_id or not m_name:
                continue
            gid = m_id.group(1)
            gnm = m_name.group(1)
            gene_id_to_name[gid] = gnm
            gene_name_to_ids.setdefault(gnm.lower(), set()).add(gid)
    return gene_name_to_ids, gene_id_to_name

def load_hgnc_synonyms(path: str):
    df = pd.read_csv(path, sep='\t', dtype=str, low_memory=False)
    cols = {c.lower(): c for c in df.columns}
    sym = cols.get('symbol') or 'symbol'
    alias = cols.get('alias_symbol')
    prev = cols.get('prev_symbol')
    mapping = {}
    def upd(tokens, canonical):
        if not isinstance(tokens, str) or not tokens.strip():
            return
        for t in re.split(r'[,|;]\s*', tokens.strip()):
            if t:
                mapping[t.lower()] = canonical
    for _, row in df.iterrows():
        canonical = str(row.get(sym, '')).strip()
        if canonical:
            mapping[canonical.lower()] = canonical
        if alias and alias in row and isinstance(row[alias], str):
            upd(row[alias], canonical)
        if prev and prev in row and isinstance(row[prev], str):
            upd(row[prev], canonical)
    return mapping