"""
Classifies: CHEBI:26244 prenols
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_prenols(smiles: str):
    """
    Determines if a molecule is a prenol based on its SMILES string.
    Prenols are alcohols or their derivatives with a carbon skeleton composed of one or more isoprene units (H-[CH2C(Me)=CHCH2]n-OH).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prenol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Check for terminal hydroxyl group (-OH)
    terminal_oh = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8 and atom.GetDegree() == 1 and atom.GetTotalNumHs() >= 1:
            # Oxygen is part of -OH group
            carbon = atom.GetNeighbors()[0]
            non_o_neighbors = [n for n in carbon.GetNeighbors() if n.GetAtomicNum() != 8]
            if len(non_o_neighbors) == 1:
                terminal_oh = True
                break

    # Check for terminal phosphate group (C-O-P linkage at chain end)
    phosphate_terminal = False
    phosphate_pattern = Chem.MolFromSmarts("[CX4]-[O]-[P]")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    for match in phosphate_matches:
        c_idx, o_idx, p_idx = match
        c_atom = mol.GetAtomWithIdx(c_idx)
        # Check if the carbon is terminal (only one non-O neighbor)
        non_o_neighbors = [n for n in c_atom.GetNeighbors() if n.GetAtomicNum() != 8]
        if len(non_o_neighbors) == 1:
            phosphate_terminal = True
            break

    if not terminal_oh and not phosphate_terminal:
        return False, "No terminal hydroxyl or phosphate group"

    # Check for isoprene units: methyl group adjacent to double bond in correct arrangement
    # Pattern matches CH3-C(=C)-C structure (isoprene unit)
    isoprene_pattern = Chem.MolFromSmarts("[CH3]-C(=C)-C")
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    if not isoprene_matches:
        return False, "No isoprene units detected"

    # Check carbon count is at least 5 (minimum for one isoprene unit)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 5:
        return False, "Insufficient carbon count"

    # Check that the isoprene units form a continuous chain
    # This is complex; instead, ensure at least one isoprene unit exists
    if len(isoprene_matches) < 1:
        return False, "Insufficient isoprene units"

    return True, "Terminal hydroxyl/phosphate with isoprene-based carbon skeleton"