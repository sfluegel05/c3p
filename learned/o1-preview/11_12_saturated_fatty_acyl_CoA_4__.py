"""
Classifies: CHEBI:84948 11,12-saturated fatty acyl-CoA(4-)
"""
"""
Classifies: 11,12-saturated fatty acyl-CoA(4-)
"""
from rdkit import Chem

def is_11_12_saturated_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 11,12-saturated fatty acyl-CoA(4-) based on its SMILES string.
    Any fatty acyl-CoA(4-) in which the 11-12 bond of the fatty acyl group is saturated.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 11,12-saturated fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define thioester pattern to find the fatty acyl-CoA linkage
    thioester_pattern = Chem.MolFromSmarts('[#6](=O)[#16]')  # C(=O)S

    # Find matches of the thioester pattern
    matches = mol.GetSubstructMatches(thioester_pattern)
    if len(matches) == 0:
        return False, "No thioester linkage found (not a fatty acyl-CoA)"

    # Assume the first match is the fatty acyl-CoA linkage
    thioester_carbon_idx = matches[0][0]

    # Function to get the longest carbon chain starting from the carbonyl carbon
    def get_longest_path(mol, start_atom_idx):
        from collections import deque

        visited = set()
        max_path = []

        queue = deque()
        queue.append((start_atom_idx, [start_atom_idx]))

        while queue:
            current_atom_idx, path = queue.popleft()

            atom = mol.GetAtomWithIdx(current_atom_idx)
            neighbors = atom.GetNeighbors()

            terminal_atom = True
            for nbr in neighbors:
                nbr_idx = nbr.GetIdx()
                if nbr_idx in path:
                    continue
                # Proceed only along carbons
                if nbr.GetAtomicNum() != 6:
                    continue
                bond = mol.GetBondBetweenAtoms(current_atom_idx, nbr_idx)
                if bond.GetBondType() not in (Chem.rdchem.BondType.SINGLE, Chem.rdchem.BondType.DOUBLE):
                    continue
                new_path = path + [nbr_idx]
                queue.append((nbr_idx, new_path))
                terminal_atom = False

            if terminal_atom:
                # Reached a terminal carbon
                if len(path) > len(max_path):
                    max_path = path

        return max_path

    # Get the acyl chain as the longest path from the thioester carbon
    acyl_chain = get_longest_path(mol, thioester_carbon_idx)

    if len(acyl_chain) < 12:
        return False, "Fatty acyl chain is too short, less than 12 carbons"

    # Get the atoms at positions 11 and 12 (C11 and C12)
    c11_atom_idx = acyl_chain[10]  # Indexing from 0
    c12_atom_idx = acyl_chain[11]

    # Get the bond between C11 and C12
    bond = mol.GetBondBetweenAtoms(c11_atom_idx, c12_atom_idx)
    if bond is None:
        return False, "No bond between C11 and C12"

    if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
        return True, "C11-C12 bond is saturated (single bond)"
    elif bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
        return False, "C11-C12 bond is unsaturated (double bond)"
    else:
        return False, "C11-C12 bond is neither single nor double"