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

    # Look for alpha-amino acid backbone: N-C-C(=O)-O
    amino_acid_pattern = Chem.MolFromSmarts('[N;!$(N-C=O)][C;!$(C=O)]([H,R0])[C](=O)[O;H1,H0-]')
    matches = mol.GetSubstructMatches(amino_acid_pattern)
    if not matches:
        return False, "Molecule is not an alpha-amino acid"

    # Identify the alpha carbon (the second atom in the pattern)
    alpha_carbons = [match[1] for match in matches]
    alpha_carbon = alpha_carbons[0]  # Assuming one amino acid per molecule

    # Get the side chain atoms (R group)
    atom_indices = set()
    backbone_atoms = set()
    backbone_atoms.update(matches[0])

    # Use a breadth-first search to find all atoms connected to the alpha carbon except backbone atoms
    to_visit = [alpha_carbon]
    visited = set()
    while to_visit:
        atom_idx = to_visit.pop()
        if atom_idx in visited or atom_idx in backbone_atoms:
            continue
        visited.add(atom_idx)
        atom = mol.GetAtomWithIdx(atom_idx)
        for neighbor in atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx not in backbone_atoms and neighbor_idx not in visited:
                to_visit.append(neighbor_idx)
        atom_indices.add(atom_idx)

    # If there are no side chain atoms, it's glycine (non-polar)
    if not atom_indices:
        return False, "Glycine has no side chain, hence non-polar"

    # Create side chain molecule
    side_chain = Chem.PathToSubmol(mol, atom_indices)
    if side_chain is None:
        return False, "Cannot extract side chain"

    # Analyze side chain for hydrogen bond donors or acceptors
    num_hbd = rdMolDescriptors.CalcNumHBD(side_chain)
    num_hba = rdMolDescriptors.CalcNumHBA(side_chain)
    if num_hbd + num_hba > 0:
        return True, "Side chain can form hydrogen bonds"
    else:
        return False, "Side chain cannot form hydrogen bonds"