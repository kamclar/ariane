def infer_variant_type(c_notation: str, p_notation: str) -> str:
    """Infer the normalized variant type used by ENIGMA rule modules."""
    import re
    c = (c_notation or "").lower()
    p = (p_notation or "").lower()

    # 5'UTR: c.-N notation (negative position before start codon, e.g. c.-10A>G)
    # Must be checked before the generic "-" in c check to avoid "intronic" misclassification
    if re.match(r"c\.-\d", c):
        return "5utr"
    # 3'UTR: c.*N notation (after stop codon, e.g. c.*10A>G)
    if "*" in c:
        return "3utr"

    if "+" in c or "-" in c:
        match = re.search(r"[+-](\d+)", c)
        if match and int(match.group(1)) <= 2:
            return "splice_site"
        return "intronic"

    if "(" in c and "del" in c:
        return "exon_deletion"
    if "(" in c and "dup" in c:
        return "exon_duplication"
    if "fs" in p or "fster" in p:
        return "frameshift"
    if "ter" in p and "fs" not in p:
        return "nonsense"
    # Initiation codon: Met at position 1 only.
    # Must match Met1 followed by a letter (amino acid), not a digit.
    # Without this guard, Met1083, Met1121 etc. would be wrongly caught.
    if re.match(r"p\.\(?met1[a-z]", p):
        return "initiation_codon"
    if p and "=" in p:
        return "synonymous"
    if "delins" in c:
        return "inframe_delins" if p else "delins"
    if "del" in c:
        return "inframe_deletion" if p else "deletion"
    if "ins" in c:
        return "inframe_insertion" if p else "insertion"
    if "dup" in c:
        return "inframe_insertion" if p else "duplication"
    if p and "=" not in p and "?" not in p and "fs" not in p and "ter" not in p:
        return "missense"
    return "unknown"
