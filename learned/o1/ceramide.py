"""
Classifies: CHEBI:17761 ceramide
"""
"""
Classifies: ceramide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_ceramide(smiles: str):
    """
    Determines if a molecule is a ceramide based on its SMILES string.
    A ceramide is an N-acyl-sphingoid base with an amide-linked fatty acid.
    The fatty acid is typically saturated or monounsaturated with 14 to 26 carbons.
    The sphingoid base is a long-chain amino alcohol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ceramide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    try:
        # Define SMARTS patterns
        amide_smarts = '[CX3](=O)[NX3][C]'  # Amide bond
        sphingoid_smarts = '[NX3][C@@H](CO)[C@H](O)CCCCCCCCCCCCCCC'  # Sphingoid base core
        fatty_acid_smarts = '[C](-[C]){12,24}-C(=O)N'  # Fatty acid chain with 14-26 carbons

        # Convert SMARTS to molecules
        amide_pattern = Chem.MolFromSmarts(amide_smarts)
        sphingoid_pattern = Chem.MolFromSmarts(sphingoid_smarts)
        fatty_acid_pattern = Chem.MolFromSmarts(fatty_acid_smarts)

        # Check for amide bond
        if not mol.HasSubstructMatch(amide_pattern):
            return False, "No amide bond found"

        # Find amide matches
        amide_matches = mol.GetSubstructMatches(amide_pattern)

        # Loop through all amide matches to find valid ceramide structure
        for amide_match in amide_matches:
            carbonyl_c_idx = amide_match[0]
            nitrogen_idx = amide_match[2]

            carbonyl_c = mol.GetAtomWithIdx(carbonyl_c_idx)
            nitrogen = mol.GetAtomWithIdx(nitrogen_idx)

            # Check for sphingoid base attached to nitrogen
            sphingoid_match = mol.HasSubstructMatch(sphingoid_pattern, useChirality=True)
            if not sphingoid_match:
                continue  # Try next amide bond

            # Get fatty acid chain length starting from carbonyl carbon
            fatty_acid_length = get_chain_length(mol, carbonyl_c_idx, exclude_idx=nitrogen_idx)
            if not 14 <= fatty_acid_length <= 26:
                continue  # Try next amide bond

            # Get sphingoid base chain length starting from nitrogen
            sphingoid_length = get_chain_length(mol, nitrogen_idx, exclude_idx=carbonyl_c_idx)
            if sphingoid_length < 12:
                continue  # Try next amide bond

            # Check for hydroxyl groups on sphingoid base
            sphingoid_oxygen_count = count_heteroatoms(mol, nitrogen_idx, atomic_num=8, exclude_idx=carbonyl_c_idx)
            if sphingoid_oxygen_count < 1:
                continue  # Try next amide bond

            # All checks passed
            return True, "Molecule is a ceramide with appropriate chain lengths and functional groups"

        # If no valid ceramide structure found
        return False, "No valid ceramide structure found in molecule"

    except Exception as e:
        return False, f"Error during processing: {e}"

def get_chain_length(mol, start_idx, exclude_idx=None):
    """
    Calculates the length of a carbon chain starting from a given atom index.

    Args:
        mol (Chem.Mol): Molecule object
        start_idx (int): Starting atom index
        exclude_idx (int): Atom index to exclude from traversal

    Returns:
        int: Number of carbon atoms in the chain
    """
    visited = set()
    queue = [start_idx]
    carbon_count = 0

    while queue:
        idx = queue.pop()
        if idx in visited or idx == exclude_idx:
            continue
        visited.add(idx)
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() == 6:
            carbon_count += 1
            for neighbor in atom.GetNeighbors():
                n_idx = neighbor.GetIdx()
                if n_idx not in visited:
                    queue.append(n_idx)
    return carbon_count

def count_heteroatoms(mol, start_idx, atomic_num, exclude_idx=None):
    """
    Counts the number of heteroatoms (e.g., oxygen) connected to the chain starting from a given atom index.

    Args:
        mol (Chem.Mol): Molecule object
        start_idx (int): Starting atom index
        atomic_num (int): Atomic number of the heteroatom to count
        exclude_idx (int): Atom index to exclude from traversal

    Returns:
        int: Number of heteroatoms found
    """
    visited = set()
    queue = [start_idx]
    heteroatom_count = 0

    while queue:
        idx = queue.pop()
        if idx in visited or idx == exclude_idx:
            continue
        visited.add(idx)
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() == atomic_num:
            heteroatom_count += 1
        for neighbor in atom.GetNeighbors():
            n_idx = neighbor.GetIdx()
            if n_idx not in visited:
                queue.append(n_idx)
    return heteroatom_count