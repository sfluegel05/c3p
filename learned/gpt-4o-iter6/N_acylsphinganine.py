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
    if not mol:
        return False, "Invalid SMILES string"

    # Look for a broader sphinganine backbone pattern:
    # adjusted pattern to allow flexible connection and less strict stereo specification
    sphinganine_pattern = Chem.MolFromSmarts("[C][C](O)[C][C](O)NC")
    if sphinganine_pattern is None:
        return False, "Failed to build sphinganine pattern"
        
    if not mol.HasSubstructMatch(sphinganine_pattern):
        return False, "No sphinganine backbone found"
    
    # Look for N-acyl amide bond pattern (-NC(=O)-), modified to recognize broader connections
    amide_pattern = Chem.MolFromSmarts("N[C](=O)")
    if amide_pattern is None:
        return False, "Failed to build amide pattern"
        
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide bond (N-acyl group) found"

    # Check for a reasonable aliphatic tail chain, only needing long enough chains
    fatty_acid_chain_pattern = Chem.MolFromSmarts("C(=O)C(C){8,}")  # Require longer chains more generally
    if fatty_acid_chain_pattern is None:
        return False, "Failed to build fatty acid chain pattern"
        
    if not mol.HasSubstructMatch(fatty_acid_chain_pattern):
        return False, "No appropriate aliphatic carbon chain found attached to the amide"

    return True, "Molecule is an N-acylsphinganine with sphinganine backbone and suitably matched N-acyl group"