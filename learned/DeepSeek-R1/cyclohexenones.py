"""
Classifies: CHEBI:48953 cyclohexenones
"""
from rdkit import Chem
from rdkit.Chem import MolFromSmiles

def is_cyclohexenones(smiles: str):
    """
    Determines if a molecule is a cyclohexenone based on its SMILES string.
    A cyclohexenone is a six-membered alicyclic ketone having one double bond in the ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cyclohexenone, False otherwise
        str: Reason for classification
    """
    mol = MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()

    for ring in rings:
        if len(ring) != 6:
            continue  # Not six-membered

        # Check for aromatic bonds in the ring
        is_aromatic = False
        ring_atoms = set(ring)
        for bond in mol.GetBonds():
            if (bond.GetBeginAtomIdx() in ring_atoms and 
                bond.GetEndAtomIdx() in ring_atoms and 
                bond.GetIsAromatic()):
                is_aromatic = True
                break
        if is_aromatic:
            continue  # Skip aromatic rings

        # Count double bonds between carbons in the ring
        ring_double_bonds = 0
        for bond in mol.GetBonds():
            a1 = bond.GetBeginAtomIdx()
            a2 = bond.GetEndAtomIdx()
            if (a1 in ring_atoms and a2 in ring_atoms and
                bond.GetBondType() == Chem.BondType.DOUBLE and
                mol.GetAtomWithIdx(a1).GetAtomicNum() == 6 and
                mol.GetAtomWithIdx(a2).GetAtomicNum() == 6):
                ring_double_bonds += 1

        if ring_double_bonds != 1:
            continue  # Must have exactly one double bond in the ring

        # Check for ketone group (C=O) attached to the ring
        has_ketone = False
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() == 6:
                for bond in atom.GetBonds():
                    if (bond.GetBondType() == Chem.BondType.DOUBLE and
                        bond.GetOtherAtom(atom).GetAtomicNum() == 8 and
                        bond.GetOtherAtom(atom).GetIdx() not in ring_atoms):
                        has_ketone = True
                        break
                if has_ketone:
                    break

        if has_ketone:
            return True, "Six-membered alicyclic ketone with one double bond in the ring"

    return False, "No matching six-membered alicyclic ketone found"