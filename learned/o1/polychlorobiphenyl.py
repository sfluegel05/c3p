"""
Classifies: CHEBI:53156 polychlorobiphenyl
"""
from rdkit import Chem

def is_polychlorobiphenyl(smiles: str):
    """
    Determines if a molecule is a polychlorobiphenyl (PCB) based on its SMILES string.
    A polychlorobiphenyl is a biphenyl compound containing between 2 and 10 chlorine atoms attached to the two benzene rings.

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

    # Define biphenyl pattern
    biphenyl_smarts = "c1ccccc1-c2ccccc2"
    biphenyl_pattern = Chem.MolFromSmarts(biphenyl_smarts)
    biphenyl_matches = mol.GetSubstructMatches(biphenyl_pattern)
    if not biphenyl_matches:
        return False, "No biphenyl core found"

    # For each biphenyl core found
    for match in biphenyl_matches:
        biphenyl_atoms = match
        biphenyl_atom_set = set(biphenyl_atoms)

        chloro_count = 0
        for atom_idx in biphenyl_atoms:
            atom = mol.GetAtomWithIdx(atom_idx)
            # Check for chlorine atoms directly attached to the biphenyl carbons
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 17:  # Chlorine atomic number
                    if neighbor.GetIdx() not in biphenyl_atom_set:
                        chloro_count += 1

        if 2 <= chloro_count <= 10:
            return True, f"Contains biphenyl core with {chloro_count} chlorine atoms attached to the benzene rings"

    return False, "No biphenyl core with 2 to 10 chlorines found"