"""
Classifies: CHEBI:72588 semisynthetic derivative
"""
from rdkit import Chem

def is_semisynthetic_derivative(smiles: str):
    """
    Determines if a molecule is a semisynthetic derivative based on its SMILES string.
    Employs more comprehensive and specific substructure patterns commonly seen in semisynthetic derivatives.

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

    # Enhanced substructure patterns
    lactone_pattern = Chem.MolFromSmarts("C1OC(=O)C=C1")  # Extended lactone structures
    beta_lactam_pattern = Chem.MolFromSmarts("C1NC(=O)[C@H]C1")  # Extended beta-lactam for complex antibiotics
    ester_modification_pattern = Chem.MolFromSmarts("[C](=O)O[C,R]")  # Ester modifications
    ether_pattern = Chem.MolFromSmarts("CCOC")  # Common ether linkages
    
    # Modified chiral/stereochemical complexity
    stereochemistry_pattern = Chem.MolFromSmarts("[C@&C]")  # Any complex stereocenter
    
    # Substructure occurrences
    matches_lactone = mol.HasSubstructMatch(lactone_pattern)
    matches_beta_lactam = mol.HasSubstructMatch(beta_lactam_pattern)
    matches_ester = mol.HasSubstructMatch(ester_modification_pattern)
    matches_ether = mol.HasSubstructMatch(ether_pattern)
    matches_stereochemistry = mol.HasSubstructMatch(stereochemistry_pattern)

    # Combine matches using an improved heuristic
    if any([matches_lactone, matches_beta_lactam, matches_ester]) and any([matches_ether, matches_stereochemistry]):
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