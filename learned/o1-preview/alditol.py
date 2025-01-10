"""
Classifies: CHEBI:17522 alditol
"""
"""
Classifies: CHEBI:18385 alditol
"""

from rdkit import Chem

def is_alditol(smiles: str):
    """
    Determines if a molecule is an alditol based on its SMILES string.
    An alditol is an acyclic polyol with the general formula HOCH2[CH(OH)]nCH2OH,
    formally derivable from an aldose by reduction of the carbonyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alditol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if molecule is acyclic (no rings)
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings"

    # Identify all carbon chains in the molecule
    chains = []
    for bond in mol.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        if atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 6:
            chains.append((atom1.GetIdx(), atom2.GetIdx()))

    # Find all paths (carbon chains) in the molecule
    from rdkit.Chem import rdmolops
    paths = rdmolops.FindAllPathsOfLengthN(mol, minPath=2, maxPath=mol.GetNumAtoms(), useBonds=True)

    # Filter paths that are continuous carbon chains
    carbon_chains = []
    for path in paths:
        path_atoms = [mol.GetAtomWithIdx(idx) for idx in path]
        if all(atom.GetAtomicNum() == 6 for atom in path_atoms):
            carbon_chains.append(path)

    # Check for alditol structure in the carbon chains
    for chain in carbon_chains:
        # Check that the chain is linear (no branches)
        if any(mol.GetAtomWithIdx(idx).GetDegree() > 3 for idx in chain):
            continue  # Skip branched chains
        
        # Check that terminal carbons are CH2OH groups
        start_atom = mol.GetAtomWithIdx(chain[0])
        end_atom = mol.GetAtomWithIdx(chain[-1])

        if not is_CH2OH(start_atom, mol):
            continue
        if not is_CH2OH(end_atom, mol):
            continue

        # Check that internal carbons are CH(OH)
        internal_atoms = [mol.GetAtomWithIdx(idx) for idx in chain[1:-1]]
        all_internal_CHOH = all(is_CHOH(atom, mol) for atom in internal_atoms)
        if not all_internal_CHOH:
            continue

        # Check for extra functional groups
        if has_extra_functional_groups(mol, chain):
            continue  # Skip chains with extra groups

        n_carbons = len(chain)
        return True, f"Contains an alditol chain of length {n_carbons} carbons"

    # No valid alditol chain found
    return False, "Molecule does not contain an alditol chain"

def is_CH2OH(atom, mol):
    """
    Checks if an atom is a terminal carbon in a CH2OH group.

    Args:
        atom (Chem.Atom): Atom to check
        mol (Chem.Mol): Molecule containing the atom

    Returns:
        bool: True if atom is CH2OH, False otherwise
    """
    if atom.GetAtomicNum() != 6:
        return False
    if atom.GetDegree() != 3:
        return False
    neighbors = atom.GetNeighbors()
    num_H = atom.GetTotalNumHs()
    num_OH = 0
    for neighbor in neighbors:
        if neighbor.GetAtomicNum() == 8:  # Oxygen
            if is_hydroxyl(neighbor, atom):
                num_OH += 1
    return num_H == 2 and num_OH == 1

def is_CHOH(atom, mol):
    """
    Checks if an atom is an internal carbon in a CH(OH) group.

    Args:
        atom (Chem.Atom): Atom to check
        mol (Chem.Mol): Molecule containing the atom

    Returns:
        bool: True if atom is CH(OH), False otherwise
    """
    if atom.GetAtomicNum() != 6:
        return False
    if atom.GetDegree() != 3:
        return False
    neighbors = atom.GetNeighbors()
    num_H = atom.GetTotalNumHs()
    num_OH = 0
    for neighbor in neighbors:
        if neighbor.GetAtomicNum() == 8:  # Oxygen
            if is_hydroxyl(neighbor, atom):
                num_OH += 1
    return num_H == 1 and num_OH == 1

def is_hydroxyl(oxygen_atom, carbon_atom):
    """
    Checks if an oxygen atom is part of a hydroxyl group attached to the given carbon.

    Args:
        oxygen_atom (Chem.Atom): Oxygen atom
        carbon_atom (Chem.Atom): Carbon atom

    Returns:
        bool: True if oxygen is part of hydroxyl group on carbon, False otherwise
    """
    if oxygen_atom.GetDegree() != 1:
        return False
    bond = carbon_atom.GetBondBetweenAtoms(carbon_atom.GetIdx(), oxygen_atom.GetIdx())
    if bond is None:
        return False
    return bond.GetBondType() == Chem.BondType.SINGLE

def has_extra_functional_groups(mol, chain_atom_indices):
    """
    Checks if the molecule has extra functional groups outside the given chain.

    Args:
        mol (Chem.Mol): Molecule to check
        chain_atom_indices (list of int): Atom indices of the chain

    Returns:
        bool: True if there are extra functional groups, False otherwise
    """
    chain_set = set(chain_atom_indices)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:
            continue  # Ignore carbons
        if atom.GetAtomicNum() == 1:
            continue  # Ignore hydrogens
        if atom.GetAtomicNum() == 8:
            # Allow only hydroxyl oxygens attached to chain carbons
            neighbors = atom.GetNeighbors()
            if len(neighbors) != 1:
                return True  # Oxygen connected to more than one atom
            neighbor = neighbors[0]
            if neighbor.GetIdx() not in chain_set:
                return True  # Oxygen attached outside the chain
            # Check if oxygen is part of hydroxyl group
            if not is_hydroxyl(atom, neighbor):
                return True
            continue
        # Any other atom indicates extra functional group
        return True
    return False