"""
Classifies: CHEBI:26244 prenols
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_prenols(smiles: str):
    """
    Determines if a molecule is a prenol based on its SMILES string.
    Prenols are alcohols with a carbon skeleton composed of one or more isoprene units (H-[CH2C(Me)=CHCH2]nOH).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prenol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Check for at least one hydroxyl group
    if not mol.HasSubstructMatch(Chem.MolFromSmarts("[OH]")):
        return False, "No hydroxyl group found"

    # Check for terminal hydroxyl (OH attached to chain-end carbon)
    terminal_oh = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8 and atom.GetDegree() == 1:  # Oxygen in -OH
            carbon = atom.GetNeighbors()[0]
            # Verify carbon has only one non-O neighbor (chain-end)
            non_o_neighbors = [n for n in carbon.GetNeighbors() if n.GetAtomicNum() != 8]
            if len(non_o_neighbors) == 1:
                terminal_oh = True
                break
    if not terminal_oh:
        return False, "No terminal hydroxyl group"

    # Check for isoprene pattern: methyl adjacent to double bond (CH3-C=C or C=C-CH3)
    isoprene_pattern1 = Chem.MolFromSmarts("[CH3]-[C]=[C]")
    isoprene_pattern2 = Chem.MolFromSmarts("[C]=[C](-[CH3])")
    if not mol.HasSubstructMatch(isoprene_pattern1) and not mol.HasSubstructMatch(isoprene_pattern2):
        return False, "No methyl-adjacent double bond (isoprene unit)"

    # Check for minimum chain length (at least one isoprene unit implies 4+ carbons)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 5:  # Smallest prenol (C5H8O) has 5 carbons (e.g., prenol CC(C)=CCO)
        return False, "Molecule too small to contain isoprene units"

    return True, "Terminal hydroxyl group with isoprene-based carbon skeleton"