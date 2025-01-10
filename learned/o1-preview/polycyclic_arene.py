"""
Classifies: CHEBI:33848 polycyclic arene
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polycyclic_arene(smiles: str):
    """
    Determines if a molecule is a polycyclic arene (polycyclic aromatic hydrocarbon) based on its SMILES string.
    A polycyclic arene is a hydrocarbon composed of fused aromatic rings.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polycyclic arene, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check that the molecule contains only carbon and hydrogen atoms
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in (6, 1):  # Atomic numbers for C and H
            return False, "Molecule contains atoms other than carbon and hydrogen"

    # Get ring information
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings < 2:
        return False, "Molecule is not polycyclic (less than two rings)"

    # Check that all rings are aromatic
    aromatic_ring_count = 0
    for ring_bond_indices in ring_info.BondRings():
        # Check if all bonds in the ring are aromatic
        if all(mol.GetBondWithIdx(bond_idx).GetIsAromatic() for bond_idx in ring_bond_indices):
            aromatic_ring_count += 1
        else:
            return False, "Not all rings are aromatic"

    if aromatic_ring_count != num_rings:
        return False, "Not all rings are aromatic"

    # Check if rings are fused (share at least two atoms)
    fused = False
    ring_atom_sets = [set(ring) for ring in ring_info.AtomRings()]
    for i, ring1 in enumerate(ring_atom_sets):
        for ring2 in ring_atom_sets[i+1:]:
            shared_atoms = ring1 & ring2
            if len(shared_atoms) >= 2:
                fused = True
                break
        if fused:
            break

    if not fused:
        return False, "Rings are not fused"

    # The molecule passed all checks
    return True, "Molecule is a polycyclic arene (polycyclic aromatic hydrocarbon)"