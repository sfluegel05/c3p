"""
Classifies: CHEBI:87659 dodecanoate ester
"""
"""
Classifies: CHEBI:<id> dodecanoate ester
"""
from rdkit import Chem

def is_dodecanoate_ester(smiles: str):
    """
    Determines if a molecule is a dodecanoate ester based on its SMILES string.
    A dodecanoate ester is an ester where the carboxylic acid component is lauric acid (dodecanoic acid).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dodecanoate ester, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define ester functional group SMARTS pattern
    ester_smarts = '[CX3](=O)[OX2][#6]'  # Ester group: C(=O)O-C
    ester_pattern = Chem.MolFromSmarts(ester_smarts)
    ester_matches = mol.GetSubstructMatches(ester_pattern)

    if not ester_matches:
        return False, "No ester groups found"

    # Define lauric acid acyl group (dodecanoyl) as SMARTS
    lauric_acid_smarts = 'C(=O)CCCCCCCCCCC'  # Carbonyl carbon followed by 11 carbons
    lauric_acid_mol = Chem.MolFromSmarts(lauric_acid_smarts)

    for match in ester_matches:
        carbonyl_c_idx = match[0]
        ester_o_idx = match[1]
        alkoxy_c_idx = match[2]

        # Check for bond between ester oxygen and alkoxy carbon
        bond = mol.GetBondBetweenAtoms(ester_o_idx, alkoxy_c_idx)
        if bond is None:
            continue  # Skip if no bond exists

        # Break the bond between ester oxygen and alkoxy carbon to get acyl fragment
        acyl_frag = Chem.FragmentOnBonds(mol, [bond.GetIdx()], addDummies=True)

        # Get fragments resulting from bond breakage
        frags = Chem.GetMolFrags(acyl_frag, asMols=True, sanitizeFrags=True)

        # Identify the acyl fragment containing the carbonyl carbon
        acyl_fragment = None
        for frag in frags:
            if frag.HasSubstructMatch(Chem.MolFromSmiles('C=O')):
                acyl_fragment = frag
                break

        if acyl_fragment is None:
            continue  # Cannot find acyl fragment, skip to next ester

        # Check if acyl fragment matches lauric acid
        if acyl_fragment.HasSubstructMatch(lauric_acid_mol):
            # Further validate that the acyl chain is linear and saturated
            if is_linear_saturated_chain(acyl_fragment, carbonyl_c_idx, expected_length=12):
                return True, "Contains lauric acid ester group"
            else:
                continue  # Acyl chain is not linear or not of correct length

    return False, "No lauric acid ester groups found"

def is_linear_saturated_chain(mol, start_idx, expected_length):
    """
    Checks if the molecular fragment is a linear saturated carbon chain of the expected length,
    starting from the specified atom index.

    Args:
        mol (Chem.Mol): Molecular fragment to check
        start_idx (int): Atom index to start traversal from
        expected_length (int): Expected number of carbon atoms including the starting carbon

    Returns:
        bool: True if the fragment is a linear saturated chain of the expected length, False otherwise
    """
    visited = set()
    to_visit = [(start_idx, 0)]  # (atom index, current length)
    max_length = 0

    while to_visit:
        current_idx, length = to_visit.pop()
        if current_idx in visited:
            continue
        visited.add(current_idx)

        atom = mol.GetAtomWithIdx(current_idx)
        if atom.GetAtomicNum() != 6:
            continue  # Not a carbon atom

        # Update maximum chain length found
        if length + 1 > max_length:
            max_length = length + 1

        neighbors = [nbr.GetIdx() for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if len(neighbors) > 2:
            return False  # Branching detected, not a linear chain

        for nbr_idx in neighbors:
            if nbr_idx not in visited:
                to_visit.append((nbr_idx, length + 1))

    # Check if the maximum chain length matches the expected length
    return max_length == expected_length