"""
Classifies: CHEBI:26244 prenols
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_prenols(smiles: str):
    """
    Determines if a molecule is a prenol based on its SMILES string.
    Prenols are alcohols or their derivatives with a carbon skeleton composed of one or more isoprene units (H-[CH2C(Me)=CHCH2]nOH).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prenol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Check for terminal hydroxyl or phosphate group
    terminal_oh = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8 and atom.GetDegree() == 1:  # Oxygen in -OH
            carbon = atom.GetNeighbors()[0]
            non_o_neighbors = [n for n in carbon.GetNeighbors() if n.GetAtomicNum() != 8]
            if len(non_o_neighbors) == 1:
                terminal_oh = True
                break

    # Check for phosphate group (C-O-P linkage)
    phosphate_pattern = Chem.MolFromSmarts("[C]-O-P")
    has_phosphate = mol.HasSubstructMatch(phosphate_pattern)

    if not terminal_oh and not has_phosphate:
        return False, "No terminal hydroxyl or phosphate group"

    # Improved isoprene unit check: methyl group adjacent to double bond in correct position
    isoprene_pattern = Chem.MolFromSmarts("[CH3]-[C]=C")
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    if not isoprene_matches:
        return False, "No isoprene units detected"

    # Check for minimum chain length (at least one isoprene unit implies 5+ carbons)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 5:
        return False, "Insufficient carbon count"

    # Additional check for consecutive isoprene units by counting methyl-double bond patterns
    # At least one is required, more indicates longer chains
    if len(isoprene_matches) < 1:
        return False, "Insufficient isoprene units"

    return True, "Terminal hydroxyl/phosphate with isoprene-based carbon skeleton"