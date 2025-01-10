"""
Classifies: CHEBI:53156 polychlorobiphenyl
"""
"""
Classifies: CHEBI:60912 polychlorobiphenyl
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polychlorobiphenyl(smiles: str):
    """
    Determines if a molecule is a polychlorobiphenyl based on its SMILES string.
    A polychlorobiphenyl is a biphenyl compound containing between 2 and 10 chlorine atoms
    attached to the two benzene rings.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polychlorobiphenyl, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for biphenyl structure (two benzene rings connected by a single bond)
    biphenyl_pattern = Chem.MolFromSmarts("c1ccccc1-c2ccccc2")
    biphenyl_matches = mol.GetSubstructMatches(biphenyl_pattern)
    
    if not biphenyl_matches:
        return False, "No biphenyl structure found"

    # Get the atoms in the biphenyl core
    biphenyl_atoms = set()
    for match in biphenyl_matches:
        biphenyl_atoms.update(match)

    # Count chlorine atoms attached to biphenyl rings
    chlorine_count = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 17:  # Chlorine
            # Check if chlorine is attached to a biphenyl ring atom
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() in biphenyl_atoms:
                    chlorine_count += 1
                    break

    # Check if the number of chlorine atoms is between 2 and 10
    if chlorine_count < 2 or chlorine_count > 10:
        return False, f"Found {chlorine_count} chlorine atoms attached to biphenyl rings, need between 2 and 10"

    # Additional check: ensure biphenyl is the main structure
    # Count non-biphenyl atoms
    non_biphenyl_atoms = 0
    for atom in mol.GetAtoms():
        if atom.GetIdx() not in biphenyl_atoms:
            non_biphenyl_atoms += 1

    # If more than 2 non-biphenyl atoms, likely not a simple polychlorobiphenyl
    if non_biphenyl_atoms > 2:
        return False, "Molecule has too many atoms outside biphenyl core"

    return True, f"Contains biphenyl structure with {chlorine_count} chlorine atoms attached to rings"