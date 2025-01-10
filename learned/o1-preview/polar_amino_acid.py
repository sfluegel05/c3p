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
    amino_acid_pattern = Chem.MolFromSmarts('[N;H1,H2][C;X4][C](=O)[O;H1,H0-]')
    matches = mol.GetSubstructMatches(amino_acid_pattern)
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

    # Identify side chain atoms connected to the alpha carbon
    side_chain_atoms = set()
    to_visit = [alpha_c_idx]
    visited = set(backbone_atoms)
    while to_visit:
        current_idx = to_visit.pop()
        if current_idx in visited:
            continue
        visited.add(current_idx)
        atom = mol.GetAtomWithIdx(current_idx)
        for neighbor in atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx not in visited and neighbor_idx not in backbone_atoms:
                to_visit.append(neighbor_idx)
                side_chain_atoms.add(neighbor_idx)

    # If there are no side chain atoms, it's glycine (non-polar)
    if not side_chain_atoms:
        return False, "Glycine has no side chain, hence non-polar"

    # Create side chain molecule
    side_chain_atom_list = list(side_chain_atoms)
    side_chain = Chem.PathToSubmol(mol, side_chain_atom_list)
    if side_chain is None or side_chain.GetNumAtoms() == 0:
        return False, "Cannot extract side chain"

    # Analyze side chain for hydrogen bond donors or acceptors
    num_hbd = rdMolDescriptors.CalcNumHBD(side_chain)
    num_hba = rdMolDescriptors.CalcNumHBA(side_chain)
    if num_hbd + num_hba > 0:
        return True, "Side chain can form hydrogen bonds"
    else:
        return False, "Side chain cannot form hydrogen bonds"