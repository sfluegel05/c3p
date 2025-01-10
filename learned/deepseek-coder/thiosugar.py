"""
Classifies: CHEBI:73754 thiosugar
"""
"""
Classifies: CHEBI:37671 thiosugar
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_thiosugar(smiles: str):
    """
    Determines if a molecule is a thiosugar based on its SMILES string.
    A thiosugar is a carbohydrate derivative where one or more oxygens or hydroxyl groups
    are replaced by sulfur or -SR groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a thiosugar, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule has a sugar-like structure (ring with multiple hydroxyl groups)
    sugar_pattern = Chem.MolFromSmarts("[C;H1,H2][OH]")  # Pattern for a carbon with a hydroxyl group
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if len(sugar_matches) < 2:  # At least 2 hydroxyl groups for a sugar-like structure
        return False, "Not enough hydroxyl groups for a sugar-like structure"

    # Check for sulfur atoms in the molecule
    sulfur_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16]
    if not sulfur_atoms:
        return False, "No sulfur atoms found"

    # Check if sulfur is attached to any part of the molecule, not just the sugar backbone
    sulfur_attached = False
    for sulfur in sulfur_atoms:
        neighbors = sulfur.GetNeighbors()
        for neighbor in neighbors:
            if neighbor.GetAtomicNum() == 6:  # Carbon atom
                sulfur_attached = True
                break
        if sulfur_attached:
            break

    if not sulfur_attached:
        return False, "Sulfur not attached to the molecule"

    # Check for the presence of a ring structure (carbohydrate-like)
    ring_info = mol.GetRingInfo()
    if not ring_info.AtomRings():
        return False, "No ring structure found"

    # Check if sulfur is replacing an oxygen or hydroxyl group in the sugar-like structure
    sulfur_replacement = False
    for sulfur in sulfur_atoms:
        for neighbor in sulfur.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # Carbon atom
                # Check if the carbon is part of a sugar-like structure
                if any(match for match in sugar_matches if neighbor.GetIdx() in match):
                    sulfur_replacement = True
                    break
        if sulfur_replacement:
            break

    if not sulfur_replacement:
        return False, "Sulfur is not replacing an oxygen or hydroxyl group in the sugar-like structure"

    return True, "Contains a sugar-like structure with sulfur replacing an oxygen or hydroxyl group"