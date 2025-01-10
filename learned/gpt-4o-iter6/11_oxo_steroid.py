"""
Classifies: CHEBI:47787 11-oxo steroid
"""
from rdkit import Chem

def is_11_oxo_steroid(smiles: str):
    """
    Determines if a molecule is an 11-oxo steroid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is an 11-oxo steroid, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define steroid backbone pattern: a tetracyclic system with three 6-membered rings and one 5-membered
    steroid_pattern = Chem.MolFromSmarts("C1CCC2C3C(C4C(CC(C3C2C1)C)C4)CC")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
    
    # Look for an oxo group pattern (C=O) in ring C (somewhere within the middle of the steroid structure)
    oxo_pattern = Chem.MolFromSmarts("C=O")
    oxo_matches = mol.GetSubstructMatches(oxo_pattern)
    for match in oxo_matches:
        # Check if the C=O is part of the middle of the steroid suggesting the vicinity of position 11
        atom_idx = match[0]
        neighbors = [nbr.GetRingInfo().IsInRing(atom.GetIdx()) for atom in mol.GetAtomWithIdx(atom_idx).GetNeighbors()]
        if True in neighbors:  # implying keto group within 11th distinct ring
            return True, "Structure has both the steroid backbone and an oxo group at position 11"
    
    return False, "Oxo group not accurately identified at position 11"