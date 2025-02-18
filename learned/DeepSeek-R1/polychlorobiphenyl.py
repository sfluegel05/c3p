"""
Classifies: CHEBI:53156 polychlorobiphenyl
"""
"""
Classifies: polychlorobiphenyl (CHEBI:xxxxx)
"""
from rdkit import Chem

def is_polychlorobiphenyl(smiles: str):
    """
    Determines if a molecule is a polychlorobiphenyl based on its SMILES string.
    A polychlorobiphenyl is a biphenyl compound with 2-10 chlorine atoms attached to the benzene rings.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polychlorobiphenyl, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Check for biphenyl structure (two benzene rings connected by a single bond)
    biphenyl_pattern = Chem.MolFromSmarts("[c]1[c][c][c][c][c]1-[c]2[c][c][c][c][c]2")
    matches = mol.GetSubstructMatches(biphenyl_pattern)
    if not matches:
        return False, "No biphenyl structure found"
    
    # Collect all atoms in the biphenyl substructure (using first match)
    biphenyl_atoms = set(matches[0])
    
    # Count chlorine atoms attached to biphenyl rings
    cl_count = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'Cl':
            # Check if any neighbor is part of the biphenyl structure
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() in biphenyl_atoms:
                    cl_count += 1
                    break  # Count each Cl once
    
    if 2 <= cl_count <= 10:
        return True, f"Polychlorobiphenyl with {cl_count} chlorine atoms"
    else:
        return False, f"Biphenyl found but chlorine count ({cl_count}) is not 2-10"