import numpy as np
import networkx as nx
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, rdFingerprintGenerator
from sklearn.metrics import precision_score, recall_score, f1_score
from itertools import product
from random import shuffle


def get_morgan_fingerprint(smiles, n_bits=2048):
    """
    Generates a Morgan fingerprint (a type of molecular fingerprint) for a given SMILES string.

    Args:
        smiles (str): The SMILES string of the molecule.
        n_bits (int): The number of bits in the fingerprint.

    Returns:
        The molecular fingerprint.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=n_bits)
    return mfpgen.GetFingerprint(mol)

def calculate_tanimoto_similarity(fp1, fp2):
    """
    Calculates the Tanimoto similarity coefficient between two molecular fingerprints.

    Args:
        fp1: Fingerprint 1.
        fp2: Fingerprint 2.

    Returns:
        float: The Tanimoto similarity score.
    """
    if fp1 is None or fp2 is None:
        return 0.0
    return DataStructs.TanimotoSimilarity(fp1, fp2)
    #return nx.jaccard_coefficient(fp1, fp2)  # RDKit Tanimoto is equivalent to Jaccard for binary vectors


def build_dti_network(interactions, drugs, targets):
    """
    Builds a bipartite graph representing drug-target interactions using NetworkX.

    Args:
        interactions (list of tuples): A list of (drug, target) interaction pairs.
        drugs (list): A list of all unique drug identifiers.
        targets (list): A list of all unique target identifiers.

    Returns:
        networkx.Graph: A bipartite graph with drugs and targets as nodes.
    """
    # Create an empty graph
    G = nx.Graph()

    # Add nodes for drugs and targets with a 'type' attribute to mark their set
    G.add_nodes_from(drugs, bipartite=0, type='drug')
    G.add_nodes_from(targets, bipartite=1, type='target')

    # Add edges for each interaction
    G.add_edges_from(interactions)

    return G


def predict_interaction(drug_smiles, target_id, all_drug_smiles, known_interactions, k=3):
    """
    Performs a simple neighborhood-based prediction for a single drug-target pair.

    Args:
        drug_smiles (str): The SMILES string of the drug to predict for.
        target_id (str): The ID of the target to predict for.
        all_drug_smiles (dict): Dictionary of all drug IDs and their SMILES. 
        known_interactions (list of tuples): All known drug-target interactions.
        k (int): The number of most similar neighbors to consider.

    Returns:
        int: 1 if an interaction is predicted, 0 otherwise.
    """
    drug_to_predict_id = list(all_drug_smiles.keys())[list(all_drug_smiles.values()).index(drug_smiles)]

    # 1. Calculate similarity between the new drug and all known drugs
    fp_to_predict = get_morgan_fingerprint(drug_smiles)
    if fp_to_predict is None:
        return 0

    similarities = {}
    for other_drug_id, other_smiles in all_drug_smiles.items():
        if other_drug_id == drug_to_predict_id:
            continue
        other_fp = get_morgan_fingerprint(other_smiles)
        similarities[other_drug_id] = calculate_tanimoto_similarity(fp_to_predict, other_fp)

    # 2. Find the 'k' most similar drugs
    most_similar_drugs = sorted(similarities, key=similarities.get, reverse=True)[:k]  

    # 3. Check if any of these similar drugs are known to interact with the target
    for similar_drug_id in most_similar_drugs:
        if (similar_drug_id, target_id) in known_interactions:
            # Prediction: If a similar drug interacts, predict a positive interaction.
            return 1

    # If no similar drugs interact with the target, predict a negative interaction.
    return 0


def evaluate_model(predictions, true_labels):
    """
    Evaluates the model's performance using precision, recall, and F1-score.

    Args:
        predictions (list): A list of binary predictions (0 or 1).
        true_labels (list): A list of the true binary labels.
    """
    precision = precision_score(true_labels, predictions)
    recall = recall_score(true_labels, predictions)
    f1 = f1_score(true_labels, predictions)

    print(f"Precision: {precision:.4f}")
    print(f"Recall:    {recall:.4f}")
    print(f"F1-Score:  {f1:.4f}")


if __name__ == "__main__":
    # --- 1. Example Data ---
    # Simplified dataset for demonstration
    drug_smiles = {
        'DrugA': 'CC(=O)Oc1ccccc1C(=O)O',  # Aspirin
        'DrugB': 'Cc1ccc(cc1)C(C)C(=O)O',  # Ibuprofen
        'DrugC': 'Oc1ccc(cc1)C(O)=O',  # Paracetamol
        'DrugD': 'CC(C)CC1CCC(=O)NC1',  # Another molecule
        'DrugE': 'c1ccc(cc1)C(=O)O'  # Benzoic acid
    }

    targets = ['TargetX', 'TargetY', 'TargetZ']

    # Known interactions (used as "training" data)
    known_interactions = [
        ('DrugA', 'TargetX'),
        ('DrugB', 'TargetX'),
        ('DrugA', 'TargetY'),
        ('DrugD', 'TargetZ')
    ]

    # --- 2. Molecular Similarity and Network Analysis ---

    # Build the drug-target network
    print("Building drug-target network...")
    G = build_dti_network(known_interactions, list(drug_smiles.keys()), targets)

    print(f"Network created with {G.number_of_nodes()} nodes and {G.number_of_edges()} edges.")
    print(f"Nodes in the graph: {list(G.nodes)}")

    # --- 3. Neighborhood-based Prediction & Validation ---

    # Define a test set (a mix of positive and negative examples)
    # Positive example: A known interaction, but we'll 'hide' it for the prediction.
    # Negative examples: Pairs that are not known to interact.
    test_set = [
        # Positive example (we 'hide' DrugB-TargetX)
        ('DrugB', 'TargetX'),
        # Negative examples
        ('DrugA', 'TargetZ'),
        ('DrugC', 'TargetX'),
        ('DrugC', 'TargetY'),
        ('DrugC', 'TargetZ'),
        ('DrugE', 'TargetX'),
        ('DrugE', 'TargetY'),
        ('DrugE', 'TargetZ')
    ]

    # True labels for the test set
    true_labels = [1, 0, 0, 0, 0, 0, 0, 0]

    # Remove the positive test example from the known interactions to simulate 'unknown' data
    known_interactions.remove(('DrugB', 'TargetX'))

    predictions = []
    print("\nStarting neighborhood-based prediction...")
    for drug_id, target_id in test_set:
        drug_smiles_to_predict = drug_smiles[drug_id]

        # Use the prediction model (with k=3 neighbors)
        prediction = predict_interaction(
            drug_smiles_to_predict,
            target_id,
            drug_smiles,
            known_interactions,
            k=3
        )
        predictions.append(prediction)

        print(
            f"Predicting interaction for ({drug_id}, {target_id}): {'Interaction predicted' if prediction == 1 else 'No interaction predicted'}")

    print("\n--- Model Evaluation ---")
    evaluate_model(predictions, true_labels)

