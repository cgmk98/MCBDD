!pip install chembl_webresource_client pandas requests

from chembl_webresource_client.new_client import new_client
import pandas as pd

def get_drugs_with_protein_targets():
    print("Fetching approved drugs...")
    molecule = new_client.molecule
    activity = new_client.activity
    target = new_client.target

    approved_drugs = molecule.filter(max_phase=4)
    results = []

    for drug in approved_drugs:
        chembl_id = drug["molecule_chembl_id"]
        name = drug.get("pref_name") or chembl_id
        approval_date = drug.get("approval_date")
        if not approval_date:
            continue

        # Show activities
        acts = activity.filter(molecule_chembl_id=chembl_id).only(["target_chembl_id", "target_type"])
        for a in acts:
            if a.get("target_type") != "SINGLE PROTEIN":
                continue

            target_info = target.get(a["target_chembl_id"])
            if not target_info:
                continue

            # UniProt-IDs 
            for comp in target_info.get("target_components", []):
                for xref in comp.get("target_component_xrefs", []):
                    if xref.get("xref_src_db") == "UniProt":
                        uniprot_id = xref["xref_id"]
                        results.append({
                            "chembl_id": chembl_id,
                            "drug_name": name,
                            "uniprot_id": uniprot_id
                        })

    return pd.DataFrame(results)

# Run it
if __name__ == "__main__":
    df = get_drugs_with_protein_targets()
    print(f"✅ Found {len(df)} drug–protein pairs.")
    print("\n--- Sample Output ---")
    print(df.head(10))  # Show the first 10 rows

