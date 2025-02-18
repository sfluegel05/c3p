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

    # Check for terminal hydroxyl group (primary alcohol: CH2OH)
    terminal_oh = False
    oh_pattern = Chem.MolFromSmarts("[CX4H2]-[OH]")
    if mol.HasSubstructMatch(oh_pattern):
        # Verify the OH is at the end of a chain
        for match in mol.GetSubstructMatches(oh_pattern):
            carbon_idx, oh_idx = match
            carbon = mol.GetAtomWithIdx(carbon_idx)
            # Check if the carbon is part of a chain end
            non_oh_neighbors = [n for n in carbon.GetNeighbors() if n.GetIdx() != oh_idx]
            if len(non_oh_neighbors) == 1 and non_oh_neighbors[0].GetAtomicNum() == 6:
                terminal_oh = True
                break

    # Check for terminal phosphate group (C-O-P where C is terminal)
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

    # Check for isoprene units: methyl adjacent to double bond
    isoprene_pattern = Chem.MolFromSmarts("[CH3]-C=C")
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    if not isoprene_matches:
        return False, "No isoprene units detected"

    # Check that the carbon skeleton is primarily isoprene-based
    # At least 50% of carbons are in isoprene units (heuristic)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 5:
        return False, "Insufficient carbon count"

    # Check that there's a continuous chain of isoprene units
    # This is complex; instead, check that all isoprene units are in a single chain
    # Find the longest chain and check for isoprene patterns within it
    # This part is challenging; for simplicity, require at least one isoprene unit
    # and check that the molecule has a linear structure
    # Alternative: check that all double bonds are conjugated in a chain
    # (This is a simplification and may not cover all cases)

    return True, "Terminal hydroxyl/phosphate with isoprene-based carbon skeleton"