
from modules.extract_signature import extract_signature
from modules.fuse_perspectives import fuse_perspectives
from modules.generate_binders import generate_binders
from modules.run_simulations import run_simulations
from modules.validate_ethics import validate_ethics
from modules.personalize_binders import personalize_binders
from modules.exporter import export_designs

def codette_pipeline(target_input):
    # Stage 1: Extract Signature
    sig = extract_signature(target_input)

    # Stage 2: Perspective Fusion
    context = fuse_perspectives(sig)

    # Stage 3: Candidate Generation
    candidates = generate_binders(context)

    # Stage 4: Simulations
    scored = run_simulations(candidates)

    # Stage 5: Ethics Filter
    ethics_checked = validate_ethics(scored)

    # Stage 6: Personalization
    personalized = personalize_binders(ethics_checked, patient_data={
        "immune_profile": ["A*24:02", "B*27:05"],
        "metabolic_rate": 1.2,
        "prior_exposure": ["SARS-CoV-2", "Influenza-B"],
        "ancestry_profile": ["Native", "Irish"]
    })

    # Stage 7: Export
    result = export_designs(personalized)
    return result

if __name__ == "__main__":
    # Example input
    test_seq = "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFD"
    output = codette_pipeline(test_seq)
    print(output)
