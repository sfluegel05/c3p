"""
Classifies: CHEBI:62499 methyl-branched fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdqueries

def is_methyl_branched_fatty_acid(smiles: str):
    """
    Determines if a molecule is a methyl-branched fatty acid based on its SMILES string.
    A methyl-branched fatty acid is a fatty acid containing methyl branches only.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a methyl-branched fatty acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify carboxylic acid groups (COOH)
    carboxylic_acid_pattern = Chem.MolFromSmarts('C(=O)[O;H1,-1]')
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if len(carboxylic_acid_matches) != 1:
        return False, f"Found {len(carboxylic_acid_matches)} carboxylic acid groups, need exactly 1"

    # Get the carbon atom of the carboxyl group
    carboxyl_c_idx = carboxylic_acid_matches[0][0]
    carboxyl_c_atom = mol.GetAtomWithIdx(carboxyl_c_idx)

    # Traverse the main carbon chain starting from the carboxyl carbon
    main_chain_atoms = set()
    visited = set()
    stack = [nbr.GetIdx() for nbr in carboxyl_c_atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
    while stack:
        atom_idx = stack.pop()
        if atom_idx not in visited:
            visited.add(atom_idx)
            atom = mol.GetAtomWithIdx(atom_idx)
            main_chain_atoms.add(atom_idx)
            for nbr in atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                if nbr.GetAtomicNum() == 6 and nbr_idx not in visited:
                    stack.append(nbr_idx)

    # Identify branches
    is_methyl_branched = True
    for atom_idx in main_chain_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        neighbors = [nbr.GetIdx() for nbr in atom.GetNeighbors() if nbr.GetIdx() not in main_chain_atoms]
        for nbr_idx in neighbors:
            branch_size = get_branch_size(mol, nbr_idx, main_chain_atoms)
            if branch_size > 1:
                return False, f"Branch at atom {atom_idx} is larger than methyl group"
            else:
                # Check if the branch is a methyl group
                branch_atom = mol.GetAtomWithIdx(nbr_idx)
                if branch_atom.GetAtomicNum() != 6:
                    return False, f"Branch at atom {atom_idx} is not a carbon atom"
                if branch_atom.GetDegree() != 1:
                    return False, f"Branch at atom {atom_idx} is not a methyl group"

    return True, "Molecule is a methyl-branched fatty acid"

def get_branch_size(mol, start_idx, main_chain_atoms):
    """
    Calculates the size of a branch originating from a specified atom index.

    Args:
        mol (Chem.Mol): RDKit molecule object
        start_idx (int): Atom index where the branch starts
        main_chain_atoms (set): Set of atom indices in the main chain

    Returns:
        int: Number of atoms in the branch
    """
    branch_atoms = set()
    visited = set()
    stack = [start_idx]
    while stack:
        atom_idx = stack.pop()
        if atom_idx not in visited and atom_idx not in main_chain_atoms:
            visited.add(atom_idx)
            branch_atoms.add(atom_idx)
            atom = mol.GetAtomWithIdx(atom_idx)
            for nbr in atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                if nbr_idx not in visited:
                    stack.append(nbr_idx)
    return len(branch_atoms)