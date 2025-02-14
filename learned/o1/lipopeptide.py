"""
Classifies: CHEBI:46895 lipopeptide
"""
"""
Classifies: lipopeptide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_lipopeptide(smiles: str):
    """
    Determines if a molecule is a lipopeptide based on its SMILES string.
    A lipopeptide is a compound consisting of a peptide with an attached lipid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a lipopeptide, False otherwise
        str: Reason for classification
    """

    try:
        # Parse SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False, "Invalid SMILES string"

        # Identify peptide bonds (amide bonds between amino acids)
        peptide_bond_smarts = "[NX3][CX3](=O)[NX3]"  # N-C(=O)-N pattern indicative of peptide linkage
        peptide_bond_pattern = Chem.MolFromSmarts(peptide_bond_smarts)
        peptide_bonds = mol.GetSubstructMatches(peptide_bond_pattern)
        num_peptide_bonds = len(peptide_bonds)
        if num_peptide_bonds < 2:
            return False, f"Insufficient peptide bonds found ({num_peptide_bonds} found, need at least 2)"

        # Identify amino acid residues (alpha-amino acids)
        amino_acid_smarts = "[NX3][CX4H][CX3](=O)"  # N-C-C(=O) pattern
        amino_acid_pattern = Chem.MolFromSmarts(amino_acid_smarts)
        amino_acids = mol.GetSubstructMatches(amino_acid_pattern)
        num_amino_acids = len(amino_acids)
        if num_amino_acids < 2:
            return False, f"Insufficient amino acid residues found ({num_amino_acids} found, need at least 2)"

        # Identify long aliphatic chains (lipid chains)
        # We will look for aliphatic chains of 8 or more carbons
        min_chain_length = 8
        longest_chain_length = 0
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 6 and not atom.IsInRing():
                chain_length = get_aliphatic_chain_length(atom, set())
                if chain_length > longest_chain_length:
                    longest_chain_length = chain_length

        if longest_chain_length < min_chain_length:
            return False, f"No sufficient lipid chain found (longest aliphatic chain is {longest_chain_length} carbons, need at least {min_chain_length})"

        # Check connectivity between lipid chain and peptide
        # Extract indices of lipid chain atoms and peptide atoms
        lipid_atoms = get_lipid_atoms(mol, min_chain_length)
        if not lipid_atoms:
            return False, "No lipid chain matching criteria found"

        peptide_atoms = set()
        for match in peptide_bonds:
            peptide_atoms.update(match)

        # Check if lipid chain is connected to peptide
        is_connected = False
        for lipid_atom_idx in lipid_atoms:
            for peptide_atom_idx in peptide_atoms:
                path = Chem.rdmolops.GetShortestPath(mol, lipid_atom_idx, peptide_atom_idx)
                if path:  # If a path exists
                    is_connected = True
                    break
            if is_connected:
                break

        if not is_connected:
            return False, "Lipid chain is not connected to peptide moiety"

        return True, "Contains both peptide bonds and lipid chain connected appropriately, indicative of a lipopeptide"

    except Exception as e:
        return False, f"An error occurred during processing: {str(e)}"


def get_aliphatic_chain_length(atom, visited):
    """
    Recursively calculates the length of the aliphatic chain starting from the given atom.

    Args:
        atom (Chem.Atom): The starting atom
        visited (set): Set of visited atom indices to avoid cycles

    Returns:
        int: Length of the aliphatic chain
    """
    idx = atom.GetIdx()
    if idx in visited:
        return 0
    visited.add(idx)
    if atom.GetAtomicNum() != 6 or atom.IsInRing():
        return 0
    max_length = 1
    for neighbor in atom.GetNeighbors():
        if neighbor.GetAtomicNum() == 6 and not neighbor.IsInRing():
            length = 1 + get_aliphatic_chain_length(neighbor, visited.copy())
            if length > max_length:
                max_length = length
    return max_length

def get_lipid_atoms(mol, min_chain_length):
    """
    Identifies the atoms involved in lipid chains of sufficient length.

    Args:
        mol (Chem.Mol): The molecule
        min_chain_length (int): Minimum chain length to consider

    Returns:
        set: Set of atom indices involved in lipid chains
    """
    lipid_atoms = set()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and not atom.IsInRing():
            chain = []
            dfs_aliphatic_chain(atom, [], chain)
            if len(chain) >= min_chain_length:
                lipid_atoms.update([a.GetIdx() for a in chain])
    return lipid_atoms

def dfs_aliphatic_chain(atom, path, chain):
    """
    Depth-first search to find aliphatic chains.

    Args:
        atom (Chem.Atom): Current atom
        path (list): Current traversal path
        chain (list): Atoms in the current chain

    Returns:
        None
    """
    if atom in path:
        return
    if atom.GetAtomicNum() != 6 or atom.IsInRing():
        return
    path.append(atom)
    chain.append(atom)
    for neighbor in atom.GetNeighbors():
        if neighbor.GetAtomicNum() == 6 and not neighbor.IsInRing():
            dfs_aliphatic_chain(neighbor, path.copy(), chain)