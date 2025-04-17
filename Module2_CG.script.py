from chembl_webresource_client.new_client import new_client
import requests
import pandas as pd
from collections import defaultdict

# Step 1: Retrieve all approved drugs from ChEMBL and sort them
def get_approved_drugs():
    molecule = new_client.molecule
    approved_drugs = molecule.filter(max_phase=4, molecule_properties__full_mwt__isnull=False)
    
    # Convert to DataFrame
    data = []
    for d in approved_drugs:
        try:
            name = d.get('pref_name') or d['molecule_chembl_id']
            year = d.get('approval_date')
            if year:
                year = int(year[:4])
            data.append((d['molecule_chembl_id'], name, year))
        except:
            continue
    df = pd.DataFrame(data, columns=["chembl_id", "name", "approval_year"])
    df = df.dropna(subset=["approval_year"])
    df = df.sort_values(by=["approval_year", "name"])
    return df

# Step 2: Get protein targets (UniProt accession numbers) for drugs since 2019
def get_targets_for_recent_drugs(drug_df, year_cutoff=2019):
    target_dict = defaultdict(set)
    target = new_client.target
    activity = new_client.activity

    recent_drugs = drug_df[drug_df["approval_year"] >= year_cutoff]
    
    for _, row in recent_drugs.iterrows():
        chembl_id = row["chembl_id"]
        acts = activity.filter(molecule_chembl_id=chembl_id).only(["target_chembl_id", "target_type"])
        
        for a in acts:
            if a["target_type"] == "SINGLE PROTEIN":
                target_info = target.get(a["target_chembl_id"])
                if target_info and target_info["target_components"]:
                    for comp in target_info["target_components"]:
                        for xref in comp.get("target_component_xrefs", []):
                            if xref["xref_src_db"] == "UniProt":
                                uniprot_id = xref["xref_id"]
                                target_dict[chembl_id].add(uniprot_id)
    return target_dict

# Step 3: For each UniProt ID, get associated keywords using UniProt API
def get_uniprot_keywords(uniprot_ids):
    keywords_dict = {}
    base_url = "https://rest.uniprot.org/uniprotkb/"
    headers = {"Accept": "application/json"}

    for uid in uniprot_ids:
        response = requests.get(f"{base_url}{uid}", headers=headers)
        if response.status_code == 200:
            data = response.json()
            keywords = [kw["value"] for kw in data.get("keywords", [])]
            keywords_dict[uid] = keywords
        else:
            keywords_dict[uid] = []
    return keywords_dict

# Run all steps
if __name__ == "__main__":
    print("Fetching approved drugs...")
    drugs_df = get_approved_drugs()
    
    print("Fetching protein targets for drugs approved since 2019...")
    targets = get_targets_for_recent_drugs(drugs_df, year_cutoff=2019)

    all_uniprot_ids = set()
    for ids in targets.values():
        all_uniprot_ids.update(ids)

    print(f"Fetching UniProt keywords for {len(all_uniprot_ids)} proteins...")
    uniprot_keywords = get_uniprot_keywords(all_uniprot_ids)

    # Example output
    print("\n--- Sample Output ---")
    for chembl_id, proteins in list(targets.items())[:3]:
        print(f"Drug {chembl_id}:")
        for protein in proteins:
            print(f"  Protein: {protein}, Keywords: {uniprot_keywords.get(protein)}")
