import requests
import json
import time
import pandas as pd


def get_chembl_id(query_name, resource_type):
    """
    Searches the ChEMBL API for a given name and returns its ChEMBL ID.

    Args:
        query_name (str): The name of the drug or protein (e.g., "Aspirin").
        resource_type (str): The type of resource to search for, either "molecule" or "target".

    Returns:
        str: The ChEMBL ID (e.g., "CHEMBL25") or None if not found.
    """
    if resource_type not in ["molecule", "target"]:
        print("Error: resource_type must be 'molecule' or 'target'.")
        return None

    base_url = "https://www.ebi.ac.uk/chembl/api/data/"+resource_type+"/search"

    params = {'q': query_name, 'limit': 1, 'format': 'json', 'organism': 'Homo sapiens'}
    print(params)

    try:
        response = requests.get(base_url, params=params)
        response.raise_for_status()  # Raises an HTTPError for bad responses
        print(response.content)
        data = response.json()

        if resource_type == "molecule":
            results_key = "molecules"
            id_key = "molecule_chembl_id"
        else:
            results_key = "targets"
            id_key = "target_chembl_id"

        if data[results_key]:
            # Return the ChEMBL ID of the first result
            return data[results_key][0][id_key]
        else:
            print(f"No {resource_type} found for name: {query_name}")
            return None

    except requests.exceptions.RequestException as e:
        print(f"Error fetching {resource_type} data: {e}")
        return None


def get_interactions(molecule_chembl_id, target_chembl_id):
    """
    Fetches bioactivity data for a specific drug-target interaction.

    Args:
        molecule_chembl_id (str): The ChEMBL ID of the drug.
        target_chembl_id (str): The ChEMBL ID of the protein target.

    Returns:
        list: A list of dictionaries, where each dictionary represents a bioactivity record.
    """
    base_url = "https://www.ebi.ac.uk/chembl/api/data/activity"

    params = {
        'molecule_chembl_id': molecule_chembl_id,
        'target_chembl_id': target_chembl_id,
        'limit': 1000,  #Set a high limit to get all relevant results
        'format': 'json'
    }

    try:
        response = requests.get(base_url, params=params)
        response.raise_for_status()
        data = response.json()

        return data.get('activities', [])

    except requests.exceptions.RequestException as e:
        print(f"Error fetching interaction data: {e}")
        return []


def main(drug_list, protein_name):
    """
    Args:
        drug_list (list): list of drugs' names  (e.g., "Aspirin").
        protein_name (str): name of the protein target (e.g., "Cyclooxygenase-2").
    """
    all_records = []
    
    # Define drug and protein names
    print(f"Searching for ChEMBL ID for target: '{protein_name}'...")
    protein_id = get_chembl_id(protein_name, "target")
    if protein_id:
        print(f"Found ChEMBL ID for {protein_name}: {protein_id}\n")
    else:
        return

    time.sleep(1)
    
    for drug_name in drug_list:
        print(f"\nProcessing drug: {drug_name}")
        print(f"Searching for ChEMBL ID for drug: '{drug_name}'...")
        drug_id = get_chembl_id(drug_name, "molecule")
        if drug_id:
            print(f"Found ChEMBL ID for {drug_name}: {drug_id}\n")
        else:
            return

        time.sleep(1)  # Small delay to avoid overwhelming the API

    
        print(f"Fetching interactions between {drug_name} ({drug_id}) and {protein_name} ({protein_id})...")
    
        interactions = get_interactions(drug_id, protein_id)
        if interactions:
            print(f"Found {len(interactions)} interactions for {drug_name} with {protein_name}.\n")
            for record in interactions:
                    all_records.append({
                        "drug_name": drug_name,
                        "drug_chembl_id": drug_id,
                        "target_name": protein_name,
                        "target_chembl_id": protein_id,
                        "canonical_smiles": record.get("canonical_smiles"),
                        "standard_type": record.get("standard_type"),
                        "standard_value": record.get("standard_value"),
                        "standard_units": record.get("standard_units"),
                        "pchembl_value": record.get("pchembl_value"),
                        "document_doi": record.get("document_doi"),
                        "source_url": f"https://www.ebi.ac.uk/chembl/g/#search_results/documents/{record.get('document_chembl_id')}"
                    })

                # Print summary for first 5
            for i, record in enumerate(interactions[:5]):
                print(f"\n--- Interaction Record {i + 1} ---")
                print(f"Type: {record.get('assay_type')}")
                print(f"Standard Type: {record.get('standard_type')}")
                print(f"Standard Value: {record.get('standard_value')}")
                print(f"Standard Units: {record.get('standard_units')}")
                print(f"pChembl Value: {record.get('pchembl_value')}")
                print(f"Document DOI: {record.get('document_doi')}")
                print(
                        f"Source Document: https://www.ebi.ac.uk/chembl/g/#search_results/documents/{record.get('document_chembl_id')}")
        else:
            print(f"No interactions found for {drug_name} with {protein_name}.\n")
            all_records.append({
                "drug_name": drug_name,
                    "drug_chembl_id": drug_id,
                    "target_name": protein_name,
                    "target_chembl_id": protein_id,
                    "canonical_smiles": None,
                    "standard_type": None,
                    "standard_value": 0,
                    "standard_units": None,
                    "pchembl_value": None,
                    "document_doi": None,
                    "source_url": None
                })
                    
                
        time.sleep(1)
    return all_records
 
    

#RUN MAIN----------------------------------
#main()
