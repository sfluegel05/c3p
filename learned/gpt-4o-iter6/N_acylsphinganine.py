"""
Classifies: CHEBI:31488 N-acylsphinganine
"""
from rdkit import Chem

def is_N_acylsphinganine(smiles: str):
    """
    Determines if a molecule is an N-acylsphinganine based on its SMILES string.
    A ceramide consisting of sphinganine where one amino hydrogen is replaced by a fatty acyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an N-acylsphinganine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Refined sphinganine backbone pattern: [C@@H](O)[C@H](CO)NC -
    # Looking for (2S,3R)-2-amino-1,3-dihydroxy backbone
    sphinganine_pattern = Chem.MolFromSmarts("[C@@H](O)[C@H](CO)NC(=O)")
    if sphinganine_pattern is None:
        return False, "Failed to build sphinganine pattern"
        
    if not mol.HasSubstructMatch(sphinganine_pattern):
        return False, "No sphinganine backbone found or correct stereochemistry"

    # Amide linkage pattern verification
    amide_pattern = Chem.MolFromSmarts("NC(=O)C")
    if amide_pattern is None:
        return False, "Failed to build amide pattern"
        
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide (N-acyl group) linkage found"

    # Aliphatic tail: An improved match that ensures sufficient chain length (10+ carbons)
    aliphatic_tail_pattern = Chem.MolFromSmarts("C(=O)CCCCCCCC[C;R0]")
    if aliphatic_tail_pattern is None:
        return False, "Failed to build aliphatic tail pattern"
        
    if not mol.HasSubstructMatch(aliphatic_tail_pattern):
        return False, "No appropriate aliphatic carbon chain found (insufficient length)"

    return True, "Molecule is an N-acylsphinganine with correct sphinganine backbone and matched N-acyl group"