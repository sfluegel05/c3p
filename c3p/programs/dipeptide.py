"""
Classifies: CHEBI:46761 dipeptide
"""
from rdkit import Chem
from rdkit.Chem import rdChemReactions

def is_dipeptide(smiles: str):
    """
    Determines if a molecule is a dipeptide based on its SMILES string.
    A dipeptide is a molecule that contains two amino acid residues connected by peptide linkages.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dipeptide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define patterns for peptide bond (amide bond between amino acids)
    peptide_bond_smarts = '[C;$(C(=O))]-[N;$(N[^$(N=C)])]'
    peptide_bond = Chem.MolFromSmarts(peptide_bond_smarts)

    # Identify peptide bonds
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond)

    # Build a graph where nodes are amino acid residues, edges are peptide bonds
    # Initialize atom-to-residue mapping
    atom_residue_map = {}
    residue_idx = 0

    # Define a pattern for the backbone nitrogen of an amino acid residue
    backbone_nitrogen = Chem.MolFromSmarts('[N;!H0;!$(N[C]=O)]')

    # Identify the backbone nitrogens to define residues
    backbone_nitrogen_matches = mol.GetSubstructMatches(backbone_nitrogen)

    for match in backbone_nitrogen_matches:
        atom_idx = match[0]
        if atom_idx not in atom_residue_map:
            residue_atoms = Chem.GetMolFragment(mol, atom_idx)
            for atom in residue_atoms:
                atom_residue_map[atom.GetIdx()] = residue_idx
            residue_idx += 1

    # Create a dictionary to represent the graph
    residue_graph = {}
    for i in range(residue_idx):
        residue_graph[i] = set()

    # Add edges between residues connected by peptide bonds
    for match in peptide_bond_matches:
        c_idx, n_idx = match
        res1 = atom_residue_map.get(c_idx, None)
        res2 = atom_residue_map.get(n_idx, None)
        if res1 is not None and res2 is not None and res1 != res2:
            residue_graph[res1].add(res2)
            residue_graph[res2].add(res1)

    # Find connected components (peptide chains)
    visited = set()
    chains = []

    for residue in residue_graph:
        if residue not in visited:
            chain = []
            stack = [residue]
            while stack:
                current = stack.pop()
                if current not in visited:
                    visited.add(current)
                    chain.append(current)
                    stack.extend(residue_graph[current] - visited)
            chains.append(chain)

    # Determine the length of each peptide chain
    chain_lengths = [len(chain) for chain in chains]

    # Check if the longest chain has exactly two residues
    if 2 in chain_lengths and max(chain_lengths) == 2:
        return True, "Contains two amino acid residues connected by peptide bonds"
    else:
        return False, "Does not contain exactly two amino acid residues connected by peptide bonds"