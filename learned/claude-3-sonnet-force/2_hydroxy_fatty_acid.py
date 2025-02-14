"""
Classifies: CHEBI:10283 2-hydroxy fatty acid
"""
"""
Classifies: CHEBI:36285 2-hydroxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_2_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxy fatty acid based on its SMILES string.
    A 2-hydroxy fatty acid is any fatty acid with a hydroxy (-OH) functional group
    in the alpha- or 2-position, and no other hydroxy groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxylic acid group (-C(=O)O)
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[O;H,-]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"

    # Find the longest carbon chain
    longest_chain = find_longest_carbon_chain(mol)
    if not longest_chain:
        return False, "No carbon chain found"

    # Look for a single hydroxy group (-OH) attached to C2 of the longest chain
    hydroxy_pattern = Chem.MolFromSmarts("[CH2][CH]([OH])[CH2]")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern, onlyMatch=longest_chain)
    if len(hydroxy_matches) != 1:
        return False, "Hydroxy group not found exclusively at C2 position of the longest chain"

    # Check for additional hydroxy groups
    if count_hydroxy_groups(mol, hydroxy_matches[0]) > 1:
        return False, "Found additional hydroxy groups beyond the one at C2"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    if c_count < 3:
        return False, "Too few carbons for a fatty acid"
    if o_count != 3:
        return False, "Must have exactly 3 oxygens (carboxyl and hydroxy groups)"

    # Additional checks (optional)
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)

    if n_rotatable < 3:
        return False, "Insufficient rotatable bonds for a fatty acid"
    if mol_wt < 100:
        return False, "Molecular weight too low for a fatty acid"

    return True, "Contains a carbon chain with a carboxylic acid and a single hydroxy group at C2"

def find_longest_carbon_chain(mol):
    """
    Finds the longest carbon chain in a molecule.

    Args:
        mol (Chem.Mol): RDKit molecule object

    Returns:
        list or None: Atom indices of the longest carbon chain, or None if no chain is found
    """
    longest_chain = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon atom
            chain = find_chain(mol, atom.GetIdx())
            if len(chain) > len(longest_chain):
                longest_chain = chain
    return longest_chain

def find_chain(mol, start_atom_idx):
    """
    Recursive function to find a carbon chain starting from a given atom index.

    Args:
        mol (Chem.Mol): RDKit molecule object
        start_atom_idx (int): Index of the starting atom

    Returns:
        list: Atom indices of the carbon chain
    """
    chain = [start_atom_idx]
    neighbors = mol.GetAtomWithIdx(start_atom_idx).GetNeighbors()
    for neighbor in neighbors:
        if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in chain:
            chain.extend(find_chain(mol, neighbor.GetIdx()))
    return list(set(chain))

def count_hydroxy_groups(mol, ignore_atom_idx=None):
    """
    Counts the number of hydroxy groups (-OH) in a molecule, optionally ignoring a specific atom.

    Args:
        mol (Chem.Mol): RDKit molecule object
        ignore_atom_idx (int, optional): Index of the atom to ignore (e.g., the C2 hydroxy group)

    Returns:
        int: Number of hydroxy groups
    """
    hydroxy_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    if ignore_atom_idx is not None:
        hydroxy_matches = [match for match in hydroxy_matches if ignore_atom_idx not in match]
    return len(hydroxy_matches)