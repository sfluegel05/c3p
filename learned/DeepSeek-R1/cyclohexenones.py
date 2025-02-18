"""
Classifies: CHEBI:48953 cyclohexenones
"""
from rdkit import Chem
from rdkit.Chem import MolFromSmiles

def is_cyclohexenones(smiles: str):
    """
    Determines if a molecule is a cyclohexenone based on its SMILES string.
    A cyclohexenone is a six-membered alicyclic ketone with at least one double bond in the ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cyclohexenone, False otherwise
        str: Reason for classification
    """
    mol = MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get all rings in the molecule
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()

    for ring in rings:
        if len(ring) != 6:
            continue  # Skip non-six-membered rings

        # Check if the ring is non-aromatic (alicyclic)
        is_aromatic = any(mol.GetBondBetweenAtoms(i, j).GetIsAromatic() for i in ring for j in ring if i < j)
        if is_aromatic:
            continue

        # Check for a ketone (C=O) in the ring
        has_ketone = False
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() == 6:
                for bond in atom.GetBonds():
                    if bond.GetBondType() == Chem.BondType.DOUBLE:
                        neighbor = bond.GetOtherAtom(atom)
                        if neighbor.GetAtomicNum() == 8 and bond.GetBeginAtomIdx() in ring and bond.GetEndAtomIdx() in ring:
                            has_ketone = True
                            break
                if has_ketone:
                    break
        if not has_ketone:
            continue

        # Check for at least one double bond in the ring (could be the ketone or another)
        has_double = any(bond.GetBondType() == Chem.BondType.DOUBLE 
                         for bond in mol.GetBonds() 
                         if bond.GetBeginAtomIdx() in ring and bond.GetEndAtomIdx() in ring)

        if has_double:
            return True, "Six-membered alicyclic ketone with at least one double bond in the ring"

    return False, "No six-membered alicyclic ketone with a double bond found"