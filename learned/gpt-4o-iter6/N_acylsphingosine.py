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
    
    # Flexible pattern that accounts for variations in sphingosine backbone
    sphingosine_pattern = Chem.MolFromSmarts("[C@@H](CO)[C@H](O)CC=CCCCCCCCCCC") 
    if not mol.HasSubstructMatch(sphingosine_pattern):
        return False, "No sphingosine backbone found"
    
    # Adjusted pattern for the N-acyl group allowing variable acyl chain length
    n_acyl_pattern = Chem.MolFromSmarts("N[C@@H](CO)C(=O)[C:1]")  # Allow for branching and variation
    if not mol.HasSubstructMatch(n_acyl_pattern):
        return False, "No N-acyl group found"

    return True, "Contains sphingosine backbone with N-linked acyl group"

# Test examples
smiles_examples = [
    "CCCCCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO)[C@H](O)\\C=C\\CCCCCCCCCCCCC",  # N-tetracosanoylsphingosine
    "CCCCCCCCCCCCCCCC(O)C(=O)N[C@@H](CO)[C@H](O)\\C=C\\CCCCCCCCCCCCC",   # N-2-hydroxystearoylsphingosine
    "C(CCCCCCCCCC)CC\\C=C\\[C@@H](O)[C@@H](NC(=O)CCCCCCCCCCCCCCCCCCCCCCCCCCCCCO)CO",  # N-(omega-hydroxytriacontanoyl)sphingosine
    "C(=C/CCCCCCCCCCCCC)\\[C@@H](O)[C@@H](NC(=O)CCCCCCCCCCCCCCCCCCC)CO",  # Example for branched
]

for smiles in smiles_examples:
    result, reason = is_N_acylsphingosine(smiles)
    print(f"SMILES: {smiles} -> {result}, Reason: {reason}")