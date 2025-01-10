"""
Classifies: CHEBI:52639 N-acylsphingosine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_N_acylsphingosine(smiles: str):
    """
    Determines if a molecule is an N-acylsphingosine based on its SMILES string.
    An N-acylsphingosine consists of a sphingosine backbone with an N-linked acyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acylsphingosine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Adjust sphingosine backbone pattern for trans double bond and correct alcohol orientation
    sphingosine_pattern = Chem.MolFromSmarts("[C@@H](CO)[C@H](O)CC=C") 
    if not mol.HasSubstructMatch(sphingosine_pattern):
        return False, "No sphingosine backbone found"
    
    # Pattern for N-acyl group: N-attached acyl group with variable chain length and connections
    n_acyl_pattern = Chem.MolFromSmarts("[NX3][C;H1,H2,H3][CX3](=O)[!!CX3H1,!!CX3R0]")
    if not mol.HasSubstructMatch(n_acyl_pattern):
        return False, "No N-acyl group found"

    return True, "Contains sphingosine backbone with N-linked acyl group"

# Example usage
smiles_examples = [
    "CCCCCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO)[C@H](O)\\C=C\\CCCCCCCCCCCCC",  # Example 1
    "CCCCCCCCCCCCCCCC(O)C(=O)N[C@@H](CO)[C@H](O)\\C=C\\CCCCCCCCCCCCC"    # Example 2
]
for smiles in smiles_examples:
    result, reason = is_N_acylsphingosine(smiles)
    print(f"SMILES: {smiles} -> {result}, Reason: {reason}")