"""
Classifies: CHEBI:76983 N-sulfonylurea
"""
from rdkit import Chem

def is_N_sulfonylurea(smiles: str):
    """
    Determines if a molecule is an N-sulfonylurea based on its SMILES string.
    N-sulfonylureas are characterized by a sulfonyl group attached to a nitrogen of a urea moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-sulfonylurea, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Specific SMARTS pattern to target N-sulfonylureas 
    # The pattern aims to identify: sulfonyl group connected to a nitrogen of urea
    nsulfonylurea_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[NX3,SX4](=O)(=O)[#6]")
    
    # Include additional variations for the pattern
    alternative_pattern1 = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[NX3][CX3](=[OX1])[SX4](=O)(=O)[#6,NX3]")
    alternative_pattern2 = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[SX4](=O)(=O)[NX3][CX3](=[OX1])[SX4](=O)(=O)#6")

    # Check if the molecule has these patterns
    match1 = mol.HasSubstructMatch(nsulfonylurea_pattern)
    match2 = mol.HasSubstructMatch(alternative_pattern1)
    match3 = mol.HasSubstructMatch(alternative_pattern2)

    if match1 or match2 or match3:
        return True, "Contains N-sulfonylurea moiety"
    else:
        return False, "Does not contain N-sulfonylurea moiety"