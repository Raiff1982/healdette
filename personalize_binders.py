
import random

# Simplified population HLA frequency references (can be expanded with real datasets)
HLA_REFERENCE = {
    "Native": ["A*24:02", "B*35:01", "C*04:01"],
    "Irish": ["A*01:01", "B*27:05", "C*07:01"]
}

# Immunological exposure impact (dummy scoring based on pathogen diversity)
def exposure_boost(sequence, exposure_list):
    hits = sum(1 for virus in exposure_list if virus.lower() in sequence.lower())
    return round(0.05 * hits, 4)

def personalize_binders(validated_input, patient_data):
    ancestry_tags = patient_data.get("ancestry_profile", ["Irish"])
    immune_profile = patient_data.get("immune_profile", [])
    exposure_history = patient_data.get("prior_exposure", [])
    metabolic_factor = float(patient_data.get("metabolic_rate", 1.0))

    personalized_output = []
    for binder in validated_input.get("validated_binders", []):
        sequence = binder["sequence"]
        base_score = (binder["stability_score"] + binder["predicted_affinity"]) / 2

        # Adjust for HLA presence
        hla_match = 0
        for tag in ancestry_tags:
            common_hlas = HLA_REFERENCE.get(tag, [])
            hla_match += sum(1 for allele in immune_profile if allele in common_hlas)

        hla_weight = 1.0 + (hla_match * 0.05)
        exposure_weight = 1.0 + exposure_boost(sequence, exposure_history)
        metabolism_weight = 1.0 / metabolic_factor  # faster metabolism = lower effective dose

        personalization_score = round(base_score * hla_weight * exposure_weight * metabolism_weight, 4)

        personalized_output.append({
            "sequence": sequence,
            "personalization_score": personalization_score,
            "ancestry_tags": ancestry_tags,
            "hla_matches": hla_match,
            "metabolic_factor": metabolic_factor,
            "exposure_weight": round(exposure_weight, 3),
            "ethics_notice": "Ancestry-aware modeling active. Logged and ethically approved by Codette's CoreConscience."
        })

    return {"personalized_binders": personalized_output}
