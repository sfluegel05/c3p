"""
Classifies: CHEBI:61778 triterpenoid saponin
"""
"""
Classifies: triterpenoid saponin
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_triterpenoid_saponin(smiles: str):
    """
    Determines if a molecule is a triterpenoid saponin based on its SMILES string.
    A triterpenoid saponin is a terpene glycoside in which the terpene moiety is a triterpenoid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triterpenoid saponin, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for triterpenoid core
    # Triterpenoids usually have 5 or more rings and are pentacyclic
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings < 5:
        return False, f"Contains {num_rings} rings, less than 5 rings required for triterpenoid core"

    # Check for fused ring system
    ssr = Chem.GetSymmSSSR(mol)
    fused_rings = 0
    for ring in ssr:
        atoms_in_ring = set(ring)
        for other_ring in ssr:
            if ring != other_ring and len(atoms_in_ring.intersection(other_ring)) > 0:
                fused_rings += 1
                break
    if fused_rings < 5:
        return False, "Insufficient fused rings for triterpenoid core"

    # Check for sugar moieties (pyranose and furanose rings)
    sugar_smarts = Chem.MolFromSmarts('*OC[C@H]1O[C@H]([C@H]([C@@H]([C@H]1O)O)O)CO')  # Glucose-like pattern
    sugars_found = len(mol.GetSubstructMatches(sugar_smarts)) > 0

    if not sugars_found:
        return False, "No sugar moieties found"

    # Check for glycosidic bond between triterpenoid and sugar(s)
    # Look for an ether linkage connecting a non-ring oxygen to the sugar ring oxygen
    glycosidic_bond_found = False
    for bond in mol.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
            if (atom1.GetAtomicNum() == 8 and atom2.GetAtomicNum() == 6) or \
               (atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 8):
                if atom1.IsInRing() and not atom2.IsInRing():
                    glycosidic_bond_found = True
                    break
                if atom2.IsInRing() and not atom1.IsInRing():
                    glycosidic_bond_found = True
                    break

    if not glycosidic_bond_found:
        return False, "No glycosidic bond connecting triterpenoid core and sugar moiety found"

    return True, "Contains triterpenoid core with sugar moiety linked via glycosidic bond"