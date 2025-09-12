
import json
import os
from datetime import datetime

def export_designs(personalized_binders, format='json', output_dir='output'):
    if format != 'json':
        raise ValueError("Only JSON format is currently supported.")

    os.makedirs(output_dir, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_file = os.path.join(output_dir, f"codette_antibody_designs_{timestamp}.json")

    with open(output_file, 'w') as f:
        json.dump(personalized_binders, f, indent=4)

    return {
        "status": "success",
        "output_file": output_file,
        "binder_count": len(personalized_binders.get("personalized_binders", []))
    }
