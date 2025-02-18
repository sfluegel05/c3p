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

    # Find potential sugar rings (5 or 6-membered with oxygen)
    sugar_rings = []
    ri = mol.GetRingInfo()
    for ring in ri.AtomRings():
        if len(ring) in [5, 6] and any(mol.GetAtomWithIdx(a).GetAtomicNum() == 8 for a in ring):
            sugar_rings.append(set(ring))

    # Find glycosidic bonds (O connecting two carbons, one from sugar ring)
    glycosidic_oxygens = set()
    for bond in mol.GetBonds():
        if bond.GetBondType() != Chem.BondType.SINGLE:
            continue
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        # Check if one atom is O and the other is C
        if {a1.GetAtomicNum(), a2.GetAtomicNum()} != {6, 8}:
            continue
        o_atom = a1 if a1.GetAtomicNum() == 8 else a2
        c_atom = a1 if a1.GetAtomicNum() == 6 else a2

        # Check if the carbon is in a sugar ring
        c_in_ring = any(c_atom.GetIdx() in r for r in sugar_rings)
        # Check if the oxygen is part of another sugar ring or connects to non-ring carbon (linear unit)
        o_in_ring = any(o_atom.GetIdx() in r for r in sugar_rings)

        if c_in_ring and not o_in_ring:
            glycosidic_oxygens.add(o_atom.GetIdx())
        elif o_in_ring:
            # Check if connected to a carbon not in the same ring
            for r in sugar_rings:
                if o_atom.GetIdx() in r and c_atom.GetIdx() not in r:
                    glycosidic_oxygens.add(o_atom.GetIdx())
                    break

    # Each glycosidic bond has one oxygen, so count unique oxygens
    glycosidic_bonds = len(glycosidic_oxygens)
    if glycosidic_bonds < 3:
        return False, f"Found {glycosidic_bonds} glycosidic bonds, need at least 3"

    # Check for four sugar units (approximated by rings + linear units with hydroxyls)
    # Split into fragments by breaking glycosidic bonds
    edmol = Chem.EditableMol(mol)
    bonds_to_remove = []
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.SINGLE:
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            if {a1.GetAtomicNum(), a2.GetAtomicNum()} == {6, 8}:
                o_atom = a1 if a1.GetAtomicNum() == 8 else a2
                if o_atom.GetIdx() in glycosidic_oxygens:
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