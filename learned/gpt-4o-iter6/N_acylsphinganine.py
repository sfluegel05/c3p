"""
Classifies: CHEBI:31488 N-acylsphinganine
"""
from rdkit import Chem

def is_N_acylsphinganine(smiles: str):
    """
    Determines if a molecule is an N-acylsphinganine based on its SMILES string,
    which is a ceramide consisting of sphinganine in which one of the amino hydrogens
    is substituted by a fatty acyl group.

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

    # Define the SMARTS pattern for sphinganine backbone:
    # 3-Hydroxy-sphinganine with proper stereochemistry (2S,3R)-2-amino-1,3-dihydroxyoctadecane
    sphinganine_pattern = Chem.MolFromSmarts("[C@@H](O)[C@H](CO)NC")
    if sphinganine_pattern is None:
        return False, "Failed to build sphinganine pattern"
        
    if not mol.HasSubstructMatch(sphinganine_pattern):
        return False, "No sphinganine backbone found"
    
    # Define the SMARTS pattern for N-acyl group via amide bond:
    amide_pattern = Chem.MolFromSmarts("N[C]=O")
    if amide_pattern is None:
        return False, "Failed to build amide pattern"
        
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide bond (N-acyl group) found"

    # Check for satisfactory aliphatic tail chain:
    # Extended carbon chain required (sufficient length not less than 10)
    aliphatic_tail_pattern = Chem.MolFromSmarts("C(=O)CCCC[C;R0]")
    if aliphatic_tail_pattern is None:
        return False, "Failed to build aliphatic tail pattern"
        
    if not mol.HasSubstructMatch(aliphatic_tail_pattern):
        return False, "No appropriate aliphatic carbon chain found"

    return True, "Molecule is an N-acylsphinganine with sphinganine backbone and suitably matched N-acyl group"