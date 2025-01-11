"""
Classifies: CHEBI:26167 polar amino acid
"""
"""
Classifies: polar amino acid
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polar_amino_acid(smiles: str):
    """
    Determines if a molecule is a polar amino acid based on its SMILES string.
    A polar amino acid is defined as any amino acid whose side chain is capable of forming 
    one or more hydrogen bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polar amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add hydrogens to accurately identify functional groups
    mol = Chem.AddHs(mol)

    # Define the amino acid backbone pattern: N-C-C(=O)-O
    backbone_pattern = Chem.MolFromSmarts('[N;H1,H2][C;X4][C](=O)[O;H1,H0-]')
    matches = mol.GetSubstructMatches(backbone_pattern)
    if not matches:
        return False, "Molecule is not an alpha-amino acid"

    # Assume the first match corresponds to the amino acid backbone
    match = matches[0]
    n_atom_idx = match[0]
    alpha_c_idx = match[1]
    c_prime_idx = match[2]
    o_atom_idx = match[3]

    # Identify backbone atoms
    backbone_atoms = set(match)

    # Also include hydrogens attached to backbone atoms
    for idx in match:
        atom = mol.GetAtomWithIdx(idx)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 1:  # Hydrogen
                backbone_atoms.add(neighbor.GetIdx())

    # Identify side chain atoms connected to the alpha carbon
    side_chain_atoms = set()
    alpha_c_atom = mol.GetAtomWithIdx(alpha_c_idx)
    for neighbor in alpha_c_atom.GetNeighbors():
        neighbor_idx = neighbor.GetIdx()
        if neighbor_idx not in backbone_atoms:
            # Start traversal from side chain atom
            atoms_to_visit = [neighbor]
            visited_atoms = set()
            while atoms_to_visit:
                current_atom = atoms_to_visit.pop()
                current_idx = current_atom.GetIdx()
                if current_idx in visited_atoms:
                    continue
                visited_atoms.add(current_idx)
                side_chain_atoms.add(current_idx)
                for next_neighbor in current_atom.GetNeighbors():
                    next_idx = next_neighbor.GetIdx()
                    if next_idx not in backbone_atoms and next_idx not in visited_atoms:
                        atoms_to_visit.append(next_neighbor)

    # If there are no side chain atoms, it's glycine (non-polar)
    if not side_chain_atoms:
        return False, "Glycine has no side chain, hence non-polar"

    # Check for hydrogen bond donors or acceptors in the side chain
    hbd_atoms = []
    hba_atoms = []

    for idx in side_chain_atoms:
        atom = mol.GetAtomWithIdx(idx)
        atomic_num = atom.GetAtomicNum()
        # Hydrogen bond donors: N-H, O-H, S-H
        if atomic_num in [7, 8, 16]:  # N, O, S
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 1:  # Hydrogen
                    hbd_atoms.append(idx)
                    break  # Only need one hydrogen to be a donor
        # Hydrogen bond acceptors: N, O, S with lone pairs
        if atomic_num in [7, 8, 16]:  # N, O, S
            # Exclude positively charged atoms (e.g., protonated amino groups)
            if atom.GetFormalCharge() >= 0:
                hba_atoms.append(idx)

    num_hbd = len(hbd_atoms)
    num_hba = len(hba_atoms)
    if num_hbd + num_hba > 0:
        return True, f"Side chain can form hydrogen bonds (HBD: {num_hbd}, HBA: {num_hba})"
    else:
        return False, "Side chain cannot form hydrogen bonds"