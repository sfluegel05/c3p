"""
Classifies: CHEBI:50126 tetrasaccharide
"""
"""
Classifies: CHEBI:166592 tetrasaccharide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Descriptors

def is_tetrasaccharide(smiles: str):
    """
    Determines if a molecule is a tetrasaccharide based on its SMILES string.
    A tetrasaccharide is an oligosaccharide comprising four monomeric monosaccharide units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrasaccharide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES"

    # Find potential glycosidic bonds (ethers connecting two carbons, at least one in a ring)
    glycosidic_bonds = 0
    sugar_rings = []
    ri = mol.GetRingInfo()
    for ring in ri.AtomRings():
        if len(ring) in [5, 6] and any(mol.GetAtomWithIdx(a).GetAtomicNum() == 8 for a in ring):
            sugar_rings.append(set(ring))

    # Find ether bonds connecting different sugar rings or to non-ring atoms (linear units)
    for bond in mol.GetBonds():
        if bond.GetBondType() != Chem.BondType.ETHER:
            continue
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        if a1.GetAtomicNum() != 8 or a2.GetAtomicNum() != 6:
            continue  # Ensure it's O connected to C
        # Check if the O is part of a glycosidic bond (connecting two sugar units)
        in_rings1 = [i for i, r in enumerate(sugar_rings) if a1.GetIdx() in r]
        in_rings2 = [i for i, r in enumerate(sugar_rings) if a2.GetIdx() in r]
        # If the O is in a ring (sugar) and connects to a non-ring C, or connects two different rings
        if (len(in_rings1) > 0 and len(in_rings2) == 0) or (len(in_rings2) > 0 and len(in_rings1) == 0) or (len(in_rings1) > 0 and len(in_rings2) > 0 and in_rings1[0] != in_rings2[0]):
            glycosidic_bonds += 1

    # Check for at least three glycosidic bonds
    if glycosidic_bonds < 3:
        return False, f"Found {glycosidic_bonds} glycosidic bonds, need at least 3"

    # Check for four sugar units (rings + linear)
    # Split into fragments by breaking glycosidic bonds
    edmol = Chem.EditableMol(mol)
    bonds_to_remove = []
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.ETHER:
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            # Check if this is a glycosidic bond as per earlier logic
            in_rings1 = [i for i, r in enumerate(sugar_rings) if a1.GetIdx() in r]
            in_rings2 = [i for i, r in enumerate(sugar_rings) if a2.GetIdx() in r]
            if (len(in_rings1) > 0 and len(in_rings2) == 0) or (len(in_rings2) > 0 and len(in_rings1) == 0) or (len(in_rings1) > 0 and len(in_rings2) > 0 and in_rings1[0] != in_rings2[0]):
                bonds_to_remove.append(bond.GetIdx())
    
    # Remove bonds and count fragments
    for bid in reversed(sorted(bonds_to_remove)):
        edmol.RemoveBond(bid)
    frags = Chem.GetMolFrags(edmol.GetMol(), asMols=True)
    if len(frags) != 4:
        return False, f"Split into {len(frags)} fragments, expected 4"

    # Check each fragment is a potential monosaccharide (has multiple OH groups)
    for frag in frags:
        oh_count = sum(1 for atom in frag.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() >= 1)
        if oh_count < 2:  # At least two hydroxyls for a monosaccharide
            return False, f"Fragment has only {oh_count} hydroxyls"

    return True, "Four sugar units connected by three glycosidic bonds"