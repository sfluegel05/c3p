"""
Classifies: CHEBI:31488 N-acylsphinganine
"""
from rdkit import Chem

def is_N_acylsphinganine(smiles: str):
    """
    Determines if a molecule is an N-acylsphinganine based on its SMILES string.
    An N-acylsphinganine is defined as a ceramide consisting of sphinganine
    in which one of the amino hydrogens is substituted by a fatty acyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acylsphinganine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # More flexible sphinganine backbone: considering possible stereochemistry and necessary groups
    sphinganine_pattern = Chem.MolFromSmarts("N[C@@H](CO)[C@@H](O)CCCCCCC")
    if not mol.HasSubstructMatch(sphinganine_pattern):
        return False, "No sphinganine backbone found"
    
    # N-acyl linkage: focus on the acyl group being at the nitrogen, allowing variation
    n_acyl_pattern = Chem.MolFromSmarts("N[C@@H](C(=O)[C,R0])")
    if not mol.HasSubstructMatch(n_acyl_pattern):
        return False, "No N-acyl substitution found"

    # Further, verify fatty acyl chain length (at least 12 carbons commonly)
    fatty_acid_pattern = Chem.MolFromSmarts("C(=O)C([R0,C]){12,}")
    if not mol.HasSubstructMatch(fatty_acid_pattern):
        return False, "No suitable long fatty acyl chain detected"

    return True, "Contains sphinganine backbone with N-acyl substitution"

# Test the function with example SMILES strings
smiles_examples = [
    "CCCCCCCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO)[C@H](O)CCCCCCCCCCCCCCC",
    "CCCCCCCCCCCCCC[C@@H](O)[C@H](CO)NC(=O)C(O)CCCCCCCCCCCC",
    "CCCCCCCCCCCCCCCCCCCCCC(O)C(=O)N[C@H](CO)[C@@H](O)CCCCCCCCCCCCCCC"
]

for smiles in smiles_examples:
    result, reason = is_N_acylsphinganine(smiles)
    print(f"SMILES: {smiles} -> {result}, {reason}")