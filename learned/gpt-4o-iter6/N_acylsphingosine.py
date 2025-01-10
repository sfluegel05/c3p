"""
Classifies: CHEBI:52639 N-acylsphingosine
"""
from rdkit import Chem

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
    
    # Improved pattern for sphingosine backbone
    sphingosine_pattern = Chem.MolFromSmarts("[C@H](O)[C@@H](N)[C@@H](CO)")  # Detailed pattern with C-N link

    # Check for presence of sphingosine backbone
    if not mol.HasSubstructMatch(sphingosine_pattern):
        return False, "No sphingosine backbone found"
    
    # Improved pattern for N-acyl group
    n_acyl_pattern = Chem.MolFromSmarts("N-C(=O)-C")  # More comprehensive acyl group check
    if not mol.HasSubstructMatch(n_acyl_pattern):
        return False, "No N-acyl group found"

    return True, "Contains sphingosine backbone with N-linked acyl group"

# Test examples
smiles_examples = [
    "CCCCCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO)[C@H](O)\\C=C\\CCCCCCCCCCCCC",  # N-tetracosanoylsphingosine
    "CCCCCCCCCCCCCCCC(O)C(=O)N[C@@H](CO)[C@H](O)\\C=C\\CCCCCCCCCCCCC",    # N-2-hydroxystearoylsphingosine
    "C(CCCCCCCCCC)CC\\C=C\\[C@@H](O)[C@@H](NC(=O)CCCCCCCCCCCCCCCCCCCCCCCCCCCCCO)CO",  # N-(omega-hydroxytriacontanoyl)sphingosine
]

for smiles in smiles_examples:
    result, reason = is_N_acylsphingosine(smiles)
    print(f"SMILES: {smiles} -> {result}, Reason: {reason}")