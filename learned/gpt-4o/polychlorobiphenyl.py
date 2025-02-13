"""
Classifies: CHEBI:53156 polychlorobiphenyl
"""
from rdkit import Chem

def is_polychlorobiphenyl(smiles: str):
    """
    Determines if a molecule is a polychlorobiphenyl based on its SMILES string.
    A polychlorobiphenyl is defined as a biphenyl compound containing between 2 and 10 chlorine atoms
    attached to the two benzene rings.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polychlorobiphenyl, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify biphenyl structure: two phenyl rings joined via a single bond
    biphenyl_pattern = Chem.MolFromSmarts("c1ccccc1-c2ccccc2")
    if not mol.HasSubstructMatch(biphenyl_pattern):
        return False, "No biphenyl structure found"

    # Count chlorine atoms attached to the biphenyl rings
    # Find atoms matching the phenyl ring pattern
    phenyl_matches = mol.GetSubstructMatches(Chem.MolFromSmarts("c1ccccc1"))
    
    if not phenyl_matches:
        return False, "No phenyl rings found"

    chlorine_count = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 17:  # Chlorine atom
            # Check connection to any phenyl ring
            if any(atom.IsInRingOfSize(6) and 
                   any(neighbor.GetIdx() in sum(phenyl_matches, ()) for neighbor in atom.GetNeighbors())
                   for phenyl in phenyl_matches):
                chlorine_count += 1

    # Check the count of chlorines
    if chlorine_count < 2 or chlorine_count > 10:
        return False, f"Contains {chlorine_count} chlorine atoms connected to biphenyl, needs between 2 and 10"

    return True, "Contains biphenyl structure with acceptable number of chlorine atoms (between 2 and 10)"