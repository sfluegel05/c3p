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
    The sphingoid base is a long-chain amino alcohol, possibly with additional hydroxyl groups or unsaturation.
    Substitutions at the primary alcohol (e.g., glycosylation) are common.

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
        amide_pattern = Chem.MolFromSmarts('C(=O)N')  # Simple amide bond
        sphingoid_base_pattern = Chem.MolFromSmarts('NCC(O)C')  # Amino alcohol
        long_chain_pattern = Chem.MolFromSmarts('CCCCCCCCCCCCCCCC')  # Long hydrocarbon chain (at least 14 carbons)

        # Check for amide bond
        amide_matches = mol.GetSubstructMatches(amide_pattern)
        if not amide_matches:
            return False, "No amide bond found"

        # Check for sphingoid base
        sphingoid_matches = mol.GetSubstructMatches(sphingoid_base_pattern)
        if not sphingoid_matches:
            return False, "No sphingoid base found"

        # Check for long hydrocarbon chain in sphingoid base
        has_long_chain = False
        for match in sphingoid_matches:
            nitrogen_idx = match[0]
            # Find connected carbon chains from nitrogen
            chain_length = get_chain_length(mol, nitrogen_idx, exclude_heteroatoms=True)
            if chain_length >= 12:
                has_long_chain = True
                break

        if not has_long_chain:
            return False, "Sphingoid base does not have a long hydrocarbon chain"

        # Check for fatty acid chain length attached to amide bond
        valid_fatty_acid = False
        for amide_match in amide_matches:
            carbonyl_c_idx = amide_match[0]
            nitrogen_idx = amide_match[2]
            # Exclude sphingoid base nitrogen
            if nitrogen_idx in [idx for match in sphingoid_matches for idx in match]:
                # Get fatty acid chain length
                chain_length = get_chain_length(mol, carbonyl_c_idx, exclude_idx=nitrogen_idx, exclude_heteroatoms=True)
                if 14 <= chain_length <= 26:
                    valid_fatty_acid = True
                    break

        if not valid_fatty_acid:
            return False, "Fatty acid chain length is not within 14 to 26 carbons"

        # All checks passed
        return True, "Molecule is classified as a ceramide"

    except Exception as e:
        return False, f"Error during processing: {e}"

def get_chain_length(mol, start_idx, exclude_idx=None, exclude_heteroatoms=False):
    """
    Calculates the length of a chain starting from a given atom index.

    Args:
        mol (Chem.Mol): Molecule object
        start_idx (int): Starting atom index
        exclude_idx (int or list): Atom index or indices to exclude from traversal
        exclude_heteroatoms (bool): Whether to exclude heteroatoms in the chain length

    Returns:
        int: Number of carbon atoms in the chain
    """
    if exclude_idx is None:
        exclude_idx = []
    elif isinstance(exclude_idx, int):
        exclude_idx = [exclude_idx]

    visited = set()
    queue = [start_idx]
    carbon_count = 0

    while queue:
        idx = queue.pop(0)
        if idx in visited or idx in exclude_idx:
            continue
        visited.add(idx)
        atom = mol.GetAtomWithIdx(idx)
        atomic_num = atom.GetAtomicNum()
        if atomic_num == 6:
            carbon_count += 1
        elif exclude_heteroatoms and atomic_num != 6:
            continue
        for neighbor in atom.GetNeighbors():
            n_idx = neighbor.GetIdx()
            if n_idx not in visited:
                queue.append(n_idx)
    return carbon_count