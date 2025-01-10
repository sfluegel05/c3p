"""
Classifies: CHEBI:72588 semisynthetic derivative
"""
from rdkit import Chem

def is_semisynthetic_derivative(smiles: str):
    """
    Determines if a molecule is a semisynthetic derivative based on its SMILES string.
    Attempts to use more relevant substructures common in semisynthetic derivatives.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule may be a semisynthetic derivative, False otherwise.
        str: Reason for classification.
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define more specific patterns indicative of semisynthetic derivatives
    lactone_pattern = Chem.MolFromSmarts("C1OC(=O)[C@@H]([C@H]1)C")  # Lactones in macrolides
    beta_lactam_pattern = Chem.MolFromSmarts("C1NC(=O)C1")  # Beta-lactam ring, common in antibiotics
    modified_chiral_pattern = Chem.MolFromSmarts("[C@&CX4]")  # Stereocenters (rechecking stereochemistry)
    unique_side_chain_pattern = Chem.MolFromSmarts("[R]C([R])C([R])[R]")  # Complex side chains

    # Substructure occurrences
    matches_lactone = mol.HasSubstructMatch(lactone_pattern)
    matches_beta_lactam = mol.HasSubstructMatch(beta_lactam_pattern)
    matches_chiral_centers = mol.HasSubstructMatch(modified_chiral_pattern)
    matches_side_chain = mol.HasSubstructMatch(unique_side_chain_pattern)

    # Combine using tuned heuristic
    if any([matches_lactone, matches_beta_lactam]) and (matches_chiral_centers or matches_side_chain):
        return True, "Significant markers indicative of semisynthetic modification are present"

    return False, "Not enough evidence for semisynthetic derivation"

# Example test case
smiles_examples = [
    "CO[C@H]1C[C@H](O[C@H]2[C@H](C)[C@@H](O[C@@H]3O[C@H](C)C[C@@H]([C@H]3OC(C)=O)N(C)C)[C@@H](C)C[C@@]3(CO3)C(=O)[C@H](C)[C@@H](OC(C)=O)[C@@H](C)[C@@H](C)OC(=O)[C@@H]2C)O[C@@H](C)[C@@H]1OC(C)=O",  # Troleandomycin
    "CCC(C)(C)C(=O)O[C@H]1C[C@@H](C)C=C2C=C[C@H](C)[C@H](CC[C@@H]3C[C@@H](O)CC(=O)O3)[C@@H]12"  # Simvastatin
]

results = [is_semisynthetic_derivative(smiles) for smiles in smiles_examples]
for idx, (result, reason) in enumerate(results):
    print(f"Example {idx + 1}: Result: {result}, Reason: {reason}")