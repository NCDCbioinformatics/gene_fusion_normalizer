import re

def normalize_fusion_string(s: str) -> str:
    if s is None:
        return ""
    s = str(s)
    s = s.replace("—","-").replace("–","-").replace("−","-")
    s = re.sub(r"\s+", "", s)
    return s

def candidate_hyphen_indices(s: str):
    return [i for i,ch in enumerate(s) if ch == "-"]

def split_fusion_smart(fusion: str, scorer):
    fusion = normalize_fusion_string(fusion)
    hyps = candidate_hyphen_indices(fusion)
    if not hyps:
        return fusion, "", -1
    best = None
    for i in hyps:
        left = fusion[:i]
        right = fusion[i+1:]
        score = scorer(left, right)
        best = max(best, (score, i, left, right)) if best else (score, i, left, right)
    _, idx, left, right = best
    return left, right, idx

def split_side_tokens(side: str, delims=",|/&"):
    if side is None or side == "":
        return []
    rx = "[" + re.escape(delims) + "]"
    toks = re.split(rx, side)
    cleaned, seen, out = [], set(), []
    for t in toks:
        t = re.sub(r"[()\[\]{}]", " ", t)
        t = re.sub(r"\s+", " ", t).strip()
        if t: cleaned.append(t)
    for t in cleaned:
        k = t.lower()
        if k not in seen:
            seen.add(k); out.append(t)
    return out