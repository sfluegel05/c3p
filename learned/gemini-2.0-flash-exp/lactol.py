"""
Classifies: CHEBI:38131 lactol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_lactol(smiles: str):
    """
    Determines if a molecule is a lactol based on its SMILES string.
    A lactol is a cyclic hemiacetal formed by intramolecular addition of a hydroxy group to an aldehydic or ketonic carbonyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lactol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    ring_info = mol.GetRingInfo()

    # 1. Check for rings
    if not any(atom.IsInRing() for atom in mol.GetAtoms()):
        return False, "Molecule does not contain a ring"

    # 2. Identify potential hemiacetal carbons
    hemiacetal_pattern = Chem.MolFromSmarts("[CX4]([OX2H1])[OX2]")
    matches = mol.GetSubstructMatches(hemiacetal_pattern)

    if not matches:
         return False, "No hemiacetal group found"

    for match in matches:
        central_carbon = mol.GetAtomWithIdx(match[0])
        hydroxyl_oxygen = mol.GetAtomWithIdx(match[1])
        ether_oxygen = mol.GetAtomWithIdx(match[2])

        # Check if all the atoms are in the same ring
        rings = ring_info.AtomRings()
        found_ring = False
        for ring in rings:
            if central_carbon.GetIdx() in ring and \
               hydroxyl_oxygen.GetIdx() in ring and \
               ether_oxygen.GetIdx() in ring:
                found_ring = True
                break
        if not found_ring:
             continue  #  hemiacetal is not completely within the same ring

        # Check number of oxygen neighbours to the central carbon, excluding hydroxyl oxygen
        other_oxygens = 0
        for neighbor in central_carbon.GetNeighbors():
            if neighbor.GetAtomicNum() == 8 and neighbor.GetIdx() != hydroxyl_oxygen.GetIdx():
                other_oxygens += 1
        if other_oxygens > 1:
            continue # Hemiacetal carbon has more than one other oxygen, not a valid lactol

        return True, "Contains a cyclic hemiacetal structure"


    return False, "No cyclic hemiacetal structure found"