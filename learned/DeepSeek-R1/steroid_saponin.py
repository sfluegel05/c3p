"""
Classifies: CHEBI:61655 steroid saponin
"""
"""
Classifies: CHEBI:??? steroid saponin (Any saponin derived from a hydroxysteroid.)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_steroid_saponin(smiles: str):
    """
    Determines if a molecule is a steroid saponin based on its SMILES string.
    A steroid saponin has a hydroxysteroid backbone with attached sugar moieties via glycosidic bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a steroid saponin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for steroid nucleus (four fused rings, cyclopentanoperhydrophenanthrene skeleton)
    # Using a SMARTS pattern for the steroid core (simplified)
    steroid_pattern = Chem.MolFromSmarts("[C@]12[C@]3[C@]4[C@](CC2)(CC[C@H]1[C@H](C)CCCC(C)C)CC[C@H]4CC3")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid nucleus detected"

    # Check for at least one hydroxyl group (-OH) attached to the steroid nucleus
    steroid_atoms = mol.GetSubstructMatch(steroid_pattern)
    hydroxyl_found = False
    for atom_idx in steroid_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() >= 1:  # OH group
                hydroxyl_found = True
                break
        if hydroxyl_found:
            break
    if not hydroxyl_found:
        return False, "No hydroxyl group on steroid nucleus"

    # Check for glycosidic bonds (O connected to steroid and to a sugar moiety)
    # Find oxygen atoms connected to the steroid nucleus and to a sugar-like structure
    glycosidic_o_found = False
    sugar_pattern = Chem.MolFromSmarts("[C][OX2][C]")  # Simplified sugar pattern (O in a ring)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8 and atom.GetDegree() == 2:  # Ether oxygen
            # Check if one side is steroid and the other is part of a sugar
            neighbors = atom.GetNeighbors()
            if len(neighbors) != 2:
                continue
            n1, n2 = neighbors
            # Check if one neighbor is in steroid nucleus and the other is in a sugar-like structure
            in_steroid = n1.GetIdx() in steroid_atoms or n2.GetIdx() in steroid_atoms
            if not in_steroid:
                continue
            # Check if the other side is part of a sugar (cyclic with multiple O's)
            sugar_mol = Chem.MolFromSmiles("C1OC(CO)C(O)C1O")  # Example: glucose
            if mol.HasSubstructMatch(sugar_mol):
                glycosidic_o_found = True
                break
            # Alternatively, check for a ring with multiple O atoms
            for n in [n1, n2]:
                if n.GetIdx() not in steroid_atoms:
                    # Check if this atom is in a ring with oxygen
                    ring_info = mol.GetRingInfo()
                    rings = ring_info.AtomRings()
                    for ring in rings:
                        if n.GetIdx() in ring:
                            o_count = sum(1 for a in ring if mol.GetAtomWithIdx(a).GetAtomicNum() == 8)
                            if o_count >= 1:  # At least one oxygen in the ring (e.g., pyranose)
                                glycosidic_o_found = True
                                break
                    if glycosidic_o_found:
                        break
            if glycosidic_o_found:
                break
    if not glycosidic_o_found:
        return False, "No glycosidic bond to sugar moiety"

    return True, "Steroid nucleus with hydroxyl group and glycosidic bond to sugar"