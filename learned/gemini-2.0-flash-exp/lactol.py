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
    if not any(ring_info.IsAtomInRing(atom.GetIdx()) for atom in mol.GetAtoms()):
        return False, "Molecule does not contain a ring"

    # 2. Identify hemiacetal carbon (C-OH and C-O-C in a ring).
    #    Look for -C(OH)-O- pattern where the C, O, and the other oxygen is in a ring.
    hemiacetal_pattern = Chem.MolFromSmarts("[CX4]([OX2H1])[OX2]1[#6]~[#6]1")
    matches = mol.GetSubstructMatches(hemiacetal_pattern)

    if not matches:
         return False, "No hemiacetal group found in a ring system"
    
    # check each match to make sure the oxygens are in the same ring and the carbon has only one other ring oxygen
    for match in matches:
        central_carbon = mol.GetAtomWithIdx(match[0])
        hydroxyl_oxygen = mol.GetAtomWithIdx(match[1])
        ether_oxygen = mol.GetAtomWithIdx(match[2])

        if not ring_info.IsAtomInRing(central_carbon.GetIdx()) or \
           not ring_info.IsAtomInRing(hydroxyl_oxygen.GetIdx()) or \
           not ring_info.IsAtomInRing(ether_oxygen.GetIdx()):
            return False, "Hemiacetal is not completely within a ring."

        
        # Check number of oxygen neighbours to the central carbon, should be exactly one other in the ring.
        other_ring_oxygens = 0
        for neighbor in central_carbon.GetNeighbors():
            if neighbor.GetAtomicNum() == 8 and ring_info.IsAtomInRing(neighbor.GetIdx()) and neighbor.GetIdx() != hydroxyl_oxygen.GetIdx():
                other_ring_oxygens += 1

        if other_ring_oxygens != 1:
            return False, "Hemiacetal carbon has more than one other ring oxygen."
    
    return True, "Contains a cyclic hemiacetal structure"