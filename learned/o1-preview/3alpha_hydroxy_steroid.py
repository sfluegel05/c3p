"""
Classifies: CHEBI:36835 3alpha-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3alpha-hydroxy steroid based on its SMILES string.
    A 3alpha-hydroxy steroid is a steroid with a hydroxyl group at position 3 in the alpha orientation.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3alpha-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for steroid nucleus (tetracyclic fused ring system)
    # Define SMARTS pattern for steroid nucleus (simplified pattern)
    steroid_smarts = '[#6]1CCC2C1CCC3C2CCC4C3CCCC4'  # Steroid core
    steroid_pattern = Chem.MolFromSmarts(steroid_smarts)
    if steroid_pattern is None:
        return None, "Invalid steroid SMARTS pattern"

    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "Steroid nucleus not found"

    # Define SMARTS pattern for 3alpha-hydroxyl group
    # Looking for a chiral carbon with S configuration attached to OH
    hydroxyl_smarts = '[C@H](O)'  # Chiral carbon with OH
    hydroxyl_pattern = Chem.MolFromSmarts(hydroxyl_smarts)
    if hydroxyl_pattern is None:
        return None, "Invalid hydroxyl SMARTS pattern"

    # Find all matches for hydroxyl pattern
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern, useChirality=True)
    if not hydroxyl_matches:
        return False, "No chiral hydroxyl group found"

    # Now check if any of these hydroxyl groups are at position 3 of the steroid nucleus
    # This is challenging because RDKit doesn't provide atom numbering according to IUPAC numbering
    # We'll approximate by checking if the hydroxylated carbon is connected to ring A of the steroid nucleus

    # Get the matched atoms for steroid nucleus
    steroid_matches = mol.GetSubstructMatch(steroid_pattern)
    if not steroid_matches:
        return False, "Steroid nucleus not properly matched"

    steroid_atoms = set(steroid_matches)

    # Identify ring A atoms (first ring in the steroid nucleus)
    ring_info = mol.GetRingInfo()
    ring_atoms = ring_info.AtomRings()
    # Find the ring that contains atoms 0 to 5 in the steroid nucleus match
    ring_A = None
    for ring in ring_atoms:
        if set(ring).issubset(steroid_atoms) and len(ring) == 6 and steroid_matches[0] in ring:
            ring_A = ring
            break
    if ring_A is None:
        return False, "Ring A of steroid nucleus not found"

    # Check if any hydroxyl group is attached to ring A
    for match in hydroxyl_matches:
        hydroxyl_carbon = match[0]
        if hydroxyl_carbon in ring_A:
            # Check if this carbon is at position 3 (approximate by degree and neighbors)
            # In steroids, position 3 carbon is connected to two other carbons in the ring
            atom = mol.GetAtomWithIdx(hydroxyl_carbon)
            if atom.GetDegree() == 3:
                # Check stereochemistry
                chiral_tag = atom.GetChiralTag()
                if chiral_tag == Chem.ChiralType.CHI_TETRAHEDRAL_CCW:
                    return True, "3alpha-hydroxy group found with correct stereochemistry"
                else:
                    return False, "Hydroxyl at position 3 but incorrect stereochemistry"
    return False, "No hydroxyl group at position 3 of steroid nucleus"