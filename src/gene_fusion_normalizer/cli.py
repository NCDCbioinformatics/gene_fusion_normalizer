import argparse, os, pandas as pd
from .version import __version__
from .utils import parse_gtf_gene_map, load_hgnc_synonyms
from .resolver import GeneResolver
from .splitter import split_fusion_smart, split_side_tokens

def read_table(path: str):
    lower = path.lower()
    if lower.endswith('.xlsx') or lower.endswith('.xls'):
        return pd.read_excel(path, engine='openpyxl')
    if lower.endswith('.csv'):
        return pd.read_csv(path)
    if lower.endswith('.tsv') or lower.endswith('.txt'):
        try: return pd.read_csv(path, sep='\t')
        except Exception: return pd.read_csv(path)
    return pd.read_csv(path, header=None, names=['fusion'])

def app():
    p = argparse.ArgumentParser(prog='gene-fusion-normalizer', description='Split fusion strings (human) and map to Ensembl gene_name & gene_id (smart rules).')
    p.add_argument('input', help='Input Excel/CSV/TSV file')
    p.add_argument('--col', default=None, help='Column containing fusion strings (default: auto-detect)')
    p.add_argument('--gtf', required=True, help='Ensembl GTF (gz/plain)')
    p.add_argument('--hgnc', help='HGNC complete set TSV for alias/prev symbols')
    p.add_argument('--split-delims', default=',|/&', help='Delimiters for multi-gene side (semicolon is preserved)')
    p.add_argument('--no-fuzzy', action='store_true', help='Disable fuzzy matching')
    p.add_argument('--explode', action='store_true', help='Explode combinations into multiple rows (default: off)')
    p.add_argument('-o','--out', default='fusion_mapped.xlsx', help='Output .xlsx or .csv')
    args = p.parse_args()

    df = read_table(args.input)
    col = args.col
    if col is None:
        for c in df.columns:
            if str(c).lower() in ('fusion','fusion_gene','original_fusion','original_gene_name','gene_fusion'):
                col = c; break
        if col is None:
            col = df.columns[0]
    else:
        target = str(col).strip().lower()
        for c in df.columns:
            if str(c).strip().lower() == target:
                col = c; break

    gene_name_to_ids, gene_id_to_name = parse_gtf_gene_map(args.gtf)
    syn = {}
    if args.hgnc and os.path.exists(args.hgnc):
        syn = load_hgnc_synonyms(args.hgnc)
    resolver = GeneResolver(gene_name_to_ids, gene_id_to_name, syn)

    def score_split(left_raw, right_raw):
        ltoks = split_side_tokens(left_raw,  delims=args.split_delims)
        rtoks = split_side_tokens(right_raw, delims=args.split_delims)
        lbest = GeneResolver.best_of([resolver.resolve_one(t, fuzzy=not args.no_fuzzy) for t in ltoks] or [resolver.resolve_one("", False)])
        rbest = GeneResolver.best_of([resolver.resolve_one(t, fuzzy=not args.no_fuzzy) for t in rtoks] or [resolver.resolve_one("", False)])
        # score: 2 if both resolve (not unmatched/special allowed), 1 if one side, 0 otherwise
        def ok(b): 
            return b is not None and b.get("match_type","unmatched") != "unmatched"
        return (1 if ok(lbest) else 0) + (1 if ok(rbest) else 0)

    rows = []
    for idx, val in enumerate(df[col].astype(str).tolist()):
        left_raw, right_raw, split_idx = split_fusion_smart(val, scorer=score_split)

        left_tokens  = split_side_tokens(left_raw,  delims=args.split_delims)
        right_tokens = split_side_tokens(right_raw, delims=args.split_delims)

        left_res  = [resolver.resolve_one(t, fuzzy=not args.no_fuzzy) for t in left_tokens] or [resolver.resolve_one("", fuzzy=False)]
        right_res = [resolver.resolve_one(t, fuzzy=not args.no_fuzzy) for t in right_tokens] or [resolver.resolve_one("", fuzzy=False)]
        left_best  = GeneResolver.best_of(left_res)
        right_best = GeneResolver.best_of(right_res)

        if args.explode and left_tokens and right_tokens:
            for l in left_res:
                for r in right_res:
                    rows.append({
                        "row_index": idx,
                        "fusion_raw": val,
                        "split_index": split_idx,
                        "left_raw": left_raw,
                        "right_raw": right_raw,
                        "left_token": l["token"],
                        "right_token": r["token"],
                        "left_symbol": l["matched_symbol"],
                        "left_ensembl_gene_id": l["ensembl_gene_id"],
                        "left_match_type": l["match_type"],
                        "right_symbol": r["matched_symbol"],
                        "right_ensembl_gene_id": r["ensembl_gene_id"],
                        "right_match_type": r["match_type"],
                        "left_tokens_all": ";".join(left_tokens),
                        "right_tokens_all": ";".join(right_tokens),
                    })
        else:
            rows.append({
                "row_index": idx,
                "fusion_raw": val,
                "split_index": split_idx,
                "left_raw": left_raw,
                "right_raw": right_raw,
                "left_token": left_best["token"] if left_best else "",
                "right_token": right_best["token"] if right_best else "",
                "left_symbol": left_best["matched_symbol"] if left_best else "",
                "left_ensembl_gene_id": left_best["ensembl_gene_id"] if left_best else "",
                "left_match_type": left_best["match_type"] if left_best else "unmatched",
                "right_symbol": right_best["matched_symbol"] if right_best else "",
                "right_ensembl_gene_id": right_best["ensembl_gene_id"] if right_best else "",
                "right_match_type": right_best["match_type"] if right_best else "unmatched",
                "left_tokens_all": ";".join(left_tokens),
                "right_tokens_all": ";".join(right_tokens),
            })

    out_df = pd.DataFrame(rows)

    if args.out.lower().endswith('.csv'):
        out_df.to_csv(args.out, index=False)
    else:
        with pd.ExcelWriter(args.out, engine='openpyxl') as xw:
            out_df.to_excel(xw, index=False, sheet_name='fusion_mapped')
    print(f"[ok] {len(out_df)} rows -> {args.out}")