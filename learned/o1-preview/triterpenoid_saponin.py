"""
Classifies: CHEBI:61778 triterpenoid saponin
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

    # Remove hydrogens for simplicity
    mol = Chem.RemoveHs(mol)

    # Check for large fused ring system characteristic of triterpenoids
    ssr = Chem.GetSymmSSSR(mol)
    ring_count = len(ssr)
    if ring_count < 4:
        return False, f"Only {ring_count} rings detected, less than 4 rings not characteristic of triterpenoids"

    # Calculate the largest ring system (number of connected rings)
    from rdkit.Chem import rdMolOps
    ri = mol.GetRingInfo()
    fused_ring_groups = []
    visited_bonds = set()
    for bond_ring in ri.BondRings():
        fused = set(bond_ring)
        bonds_to_visit = [bond_ring]
        while bonds_to_visit:
            current_ring = bonds_to_visit.pop()
            for bond_idx in current_ring:
                if bond_idx in visited_bonds:
                    continue
                visited_bonds.add(bond_idx)
                bond = mol.GetBondWithIdx(bond_idx)
                atom1 = bond.GetBeginAtom()
                atom2 = bond.GetEndAtom()
                for nbr_bond in atom1.GetBonds() + atom2.GetBonds():
                    nbr_bond_idx = nbr_bond.GetIdx()
                    if nbr_bond_idx in fused:
                        continue
                    if ri.NumBondRings(nbr_bond_idx) > 0:
                        fused.add(nbr_bond_idx)
                        bonds_to_visit.append(ri.BondRings()[nbr_bond_idx])
        fused_ring_groups.append(fused)
    # Get the largest fused ring system
    largest_fused_ring = max(fused_ring_groups, key=len)

    if len(largest_fused_ring) < 20:
        return False, "Largest fused ring system is too small to be a triterpenoid core"

    # Get atoms in the largest fused ring system
    core_atoms = set()
    for bond_idx in largest_fused_ring:
        bond = mol.GetBondWithIdx(bond_idx)
        core_atoms.add(bond.GetBeginAtomIdx())
        core_atoms.add(bond.GetEndAtomIdx())

    # Check for high carbon count in the core
    core_carbons = sum(1 for idx in core_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
    if core_carbons < 25:
        return False, f"Only {core_carbons} carbon atoms in core, less than expected for triterpenoids"

    # Detect sugar moieties (rings containing oxygen)
    sugar_rings = []
    for ring in ssr:
        atom_indices = list(ring)
        atoms = [mol.GetAtomWithIdx(idx) for idx in atom_indices]
        # Check if ring size is 5 or 6 and contains oxygen
        o_count = sum(1 for atom in atoms if atom.GetAtomicNum() == 8)
        if len(atom_indices) in [5, 6] and o_count >= 1:
            sugar_rings.append(set(atom_indices))
    if not sugar_rings:
        return False, "No sugar rings detected"

    # Check for glycosidic linkage (oxygen connecting core and sugar ring)
    core_atom_set = core_atoms
    found_linkage = False
    for sugar_ring in sugar_rings:
        shared_bonds = []
        for atom_idx in sugar_ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() == 8:
                for neighbor in atom.GetNeighbors():
                    neighbor_idx = neighbor.GetIdx()
                    if neighbor_idx in core_atom_set:
                        found_linkage = True
                        break
            if found_linkage:
                break
        if found_linkage:
            break

    if not found_linkage:
        return False, "No glycosidic linkage between core and sugar detected"

    # Check for functional groups (e.g., hydroxyl or carboxyl)
    functional_groups = ['[OX2H]', 'C(=O)[OX1H0-,OX2H1]']  # Hydroxyl and carboxyl groups
    has_functional_group = False
    for fg_smarts in functional_groups:
        fg_pattern = Chem.MolFromSmarts(fg_smarts)
        if fg_pattern and mol.HasSubstructMatch(fg_pattern):
            has_functional_group = True
            break

    if not has_functional_group:
        return False, "No hydroxyl or carboxyl groups found; unusual for triterpenoid saponins"

    return True, "Molecule is a triterpenoid saponin with a triterpenoid core and glycosidically linked sugar moiety"