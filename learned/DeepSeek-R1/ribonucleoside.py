"""
Classifies: CHEBI:18254 ribonucleoside
"""
"""
Classifies: CHEBI:46617 ribonucleoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_ribonucleoside(smiles: str):
    """
    Determines if a molecule is a ribonucleoside based on its SMILES string.
    A ribonucleoside consists of D-ribose connected to a nucleobase via a glycosidic bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ribonucleoside, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Improved D-ribose pattern with correct stereochemistry and substituent tolerance
    # Matches beta-D-ribofuranose core with C5 CH2O group
    ribose_pattern = Chem.MolFromSmarts(
        "[C@@H]1([C@H](O)[C@@H](O)[C@H](O[C@H](CO)O1)O)"
    )
    ribose_matches = mol.GetSubstructMatches(ribose_pattern)
    if not ribose_matches:
        return False, "No D-ribose sugar detected"

    # Check for phosphate groups (indicating nucleotide rather than nucleoside)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 15:  # Phosphorus
            return False, "Phosphate group present (nucleotide, not nucleoside)"

    # Check glycosidic bond - anomeric carbon connected to heteroatom (N/O)
    glycosidic_bond = False
    nucleobase_atoms = set()
    for match in ribose_matches:
        # Get anomeric carbon (C1 in the ribose pattern)
        anomeric_c = match[0]
        for neighbor in mol.GetAtomWithIdx(anomeric_c).GetNeighbors():
            if neighbor.GetAtomicNum() in [7, 8]:  # N or O
                # Track connected heteroatom and its neighbors
                glycosidic_bond = True
                nucleobase_atoms.add(neighbor.GetIdx())
                break
        if glycosidic_bond:
            break
    if not glycosidic_bond:
        return False, "No glycosidic bond to nucleobase heteroatom"

    # Check for nucleobase characteristics in connected atoms
    # Look for aromatic rings with heteroatoms (common in nucleobases)
    has_aromatic = False
    for idx in nucleobase_atoms:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetIsAromatic():
            has_aromatic = True
            break
        # Check if part of a ring with heteroatoms
        for bond in atom.GetBonds():
            neighbor = bond.GetOtherAtom(atom)
            if neighbor.GetAtomicNum() in [7, 8] and neighbor.IsInRing():
                has_aromatic = True
                break
        if has_aromatic:
            break

    if not has_aromatic:
        return False, "Connected group lacks nucleobase characteristics"

    return True, "D-ribose sugar connected to nucleobase via glycosidic bond"