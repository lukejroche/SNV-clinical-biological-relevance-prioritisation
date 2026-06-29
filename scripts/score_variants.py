import pandas as pd
import sys

inp = sys.argv[1]
out = sys.argv[2]

df = pd.read_csv(inp, sep="\t", header=None)
df.columns = [
    "chr", "pos", "ref", "alt",
    "ann", "clin_sig", "disease", "af"
]


def parse_ann(ann):
    if pd.isna(ann) or ann == ".":
        return pd.Series({
            "effect": None, "impact": None,
            "gene": None, "hgvs_c": None, "hgvs_p": None
        })
    # ANN can contain multiple transcripts — take the first (most severe)
    first = ann.split(",")[0]
    f = first.split("|")
    return pd.Series({
        "effect":  f[1]  if len(f) > 1  else None,
        "impact":  f[2]  if len(f) > 2  else None,
        "gene":    f[3]  if len(f) > 3  else None,
        "hgvs_c":  f[9]  if len(f) > 9  else None,
        "hgvs_p":  f[10] if len(f) > 10 else None,
    })

df = pd.concat([df, df["ann"].apply(parse_ann)], axis=1)


def parse_clnsig(clnsig):
    if pd.isna(clnsig) or clnsig == ".":
        return "unknown"
    # flatten any separators
    return clnsig.replace("|", "/").lower()

df["clin_sig_clean"] = df["clin_sig"].apply(parse_clnsig)

def score_af(x):
    """
    Population allele frequency from gnomAD.
    Rare variants more likely to be pathogenic (PM2 supporting evidence).
    Thresholds align with gnomAD allele frequency guidelines.
    BA1 (benign standalone) equivalent: AF > 0.05
    """
    try:
        x = float(x)
    except (ValueError, TypeError):
        # absent from gnomAD — supports PM2
        return 4, "PM2_support: absent from gnomAD"
    if x == 0:
        return 4,  "PM2_support: AF=0 in gnomAD"
    if x < 0.0001:
        return 3,  "PM2_support: AF<0.01%"
    if x < 0.001:
        return 2,  "PM2_moderate: AF<0.1%"
    if x < 0.01:
        return 0,  "common_variant: AF<1%"
    if x < 0.05:
        return -3, "BS1_support: AF>=1%"
    return -10,    "BA1_equivalent: AF>=5% likely benign"

def score_impact(impact):
    """
    SnpEff functional impact.
    HIGH impact effects (stop_gained, frameshift, splice_site) correspond
    to PVS1 evidence (loss of function) where gene mechanism supports it.
    MODERATE = missense, in-frame indel (PM4/PP3 territory).
    """
    impact_scores = {
        "HIGH":     (6,  "PVS1_candidate: high impact LOF"),
        "MODERATE": (3,  "PM4_candidate: moderate functional impact"),
        "LOW":      (0,  "low_impact: likely tolerated"),
        "MODIFIER": (-1, "modifier: non-coding or synonymous"),
    }
    return impact_scores.get(
        str(impact).upper(),
        (0, "impact_unknown")
    )

def score_effect(effect):
    """
    Specific SO effect terms from SnpEff.
    Refines the impact score for clinically important variant classes.
    Splice variants scored separately as they can be pathogenic
    even when not classified HIGH by SnpEff.
    """
    if pd.isna(effect):
        return 0, "effect_unknown"

    high_effects = {
        "stop_gained":          (6, "nonsense: premature stop"),
        "frameshift_variant":   (6, "frameshift: likely LOF"),
        "splice_acceptor_variant": (5, "splice_site: likely affects splicing"),
        "splice_donor_variant":    (5, "splice_site: likely affects splicing"),
        "start_lost":           (5, "start_lost: translation disrupted"),
        "stop_lost":            (3, "stop_lost: may extend protein"),
    }
    moderate_effects = {
        "missense_variant":         (3, "missense: requires functional assessment"),
        "inframe_insertion":        (2, "inframe_indel"),
        "inframe_deletion":         (2, "inframe_indel"),
        "splice_region_variant":    (2, "splice_region: possible splicing effect"),
        "protein_altering_variant": (2, "protein_altering"),
    }
    benign_effects = {
        "synonymous_variant":       (-1, "synonymous: likely no effect"),
        "intron_variant":           (-1, "intronic"),
        "upstream_gene_variant":    (-2, "regulatory_region"),
        "downstream_gene_variant":  (-2, "regulatory_region"),
        "intergenic_variant":       (-3, "intergenic: unlikely pathogenic"),
    }

    for effect_dict in [high_effects, moderate_effects, benign_effects]:
        if effect in effect_dict:
            return effect_dict[effect]

    return 0, "effect_not_classified"

def score_clinvar(clnsig):
    """
    ClinVar clinical significance.
    PS1: same variant previously classified pathogenic = strong evidence.
    Direct ClinVar P/LP calls are near-diagnostic in isolation.
    VUS scored neutrally — should not be used clinically without further evidence.
    """
    if "pathogenic" in clnsig and "likely" not in clnsig and "conflicting" not in clnsig:
        return 10, "PS1: ClinVar_Pathogenic"
    if "likely_pathogenic" in clnsig:
        return 7,  "PS1_support: ClinVar_LikelyPathogenic"
    if "conflicting" in clnsig:
        return 1,  "conflicting_interpretations: treat as VUS"
    if "uncertain" in clnsig or "vus" in clnsig:
        return 0,  "VUS: uncertain significance"
    if "likely_benign" in clnsig:
        return -5, "BS: ClinVar_LikelyBenign"
    if "benign" in clnsig:
        return -8, "BA1_support: ClinVar_Benign"
    return 0, "not_in_clinvar"


af_results       = df["af"].apply(score_af)
impact_results   = df["impact"].apply(score_impact)
effect_results   = df["effect"].apply(score_effect)
clinvar_results  = df["clin_sig_clean"].apply(score_clinvar)

df["af_score"],      df["af_reason"]      = zip(*af_results)
df["impact_score"],  df["impact_reason"]  = zip(*impact_results)
df["effect_score"],  df["effect_reason"]  = zip(*effect_results)
df["clinvar_score"], df["clinvar_reason"] = zip(*clinvar_results)

df["total_score"] = (
    df["af_score"] +
    df["impact_score"] +
    df["effect_score"] +
    df["clinvar_score"]
)


def assign_tier(row):
    """
    Rough tiering based on combined score and ClinVar status.
    NOT a validated ACMG classification — for prioritisation only.
    """
    if "Pathogenic" in str(row["clin_sig"]):
        return "Tier1_ClinVarConfirmed"
    if row["total_score"] >= 15:
        return "Tier1_HighPriority"
    if row["total_score"] >= 8:
        return "Tier2_ModeratePriority"
    if row["total_score"] >= 3:
        return "Tier3_LowPriority"
    return "Tier4_LikelyBenign"

df["tier"] = df.apply(assign_tier, axis=1)


output_cols = [
    "chr", "pos", "ref", "alt",
    "gene", "hgvs_c", "hgvs_p",
    "effect", "impact",
    "af", "clin_sig", "disease",
    "af_score", "impact_score", "effect_score", "clinvar_score",
    "total_score", "tier",
    "af_reason", "impact_reason", "effect_reason", "clinvar_reason"
]

df_out = df[output_cols].sort_values("total_score", ascending=False)


df_out.to_csv(out, sep="\t", index=False)

tier12 = df_out[df_out["tier"].isin(["Tier1_ClinVarConfirmed", "Tier1_HighPriority", "Tier2_ModeratePriority"])]
tier12_out = out.replace(".tsv", ".priority.tsv")
tier12.to_csv(tier12_out, sep="\t", index=False)

print(f"Total variants scored: {len(df_out)}")
print(f"Tier 1/2 variants for review: {len(tier12)}")
print(df_out["tier"].value_counts())
