"""
Classifies: CHEBI:87657 octanoate ester
"""
"""
Classifies: octanoate ester
"""

from rdkit import Chem

def is_octanoate_ester(smiles: str):
    """
    Determines if a molecule is an octanoate ester based on its SMILES string.
    An octanoate ester is any ester where the carboxylic acid component is octanoic acid (caprylic acid).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an octanoate ester, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define ester functional group pattern (Carbonyl carbon connected to ester oxygen)
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H0][#6]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)

    if not ester_matches:
        return False, "No ester groups found"

    for match in ester_matches:
        carbonyl_c_idx = match[0]  # Carbonyl carbon atom index
        ester_o_idx = match[1]     # Ester oxygen atom index

        # Traverse the acyl chain starting from the carbonyl carbon, avoiding the ester oxygen
        acyl_chain_atoms = set()
        stack = [(carbonyl_c_idx, ester_o_idx)]  # (current_atom_idx, previous_atom_idx)

        while stack:
            current_atom_idx, previous_atom_idx = stack.pop()
            if current_atom_idx in acyl_chain_atoms:
                continue  # Avoid cycles
            acyl_chain_atoms.add(current_atom_idx)
            atom = mol.GetAtomWithIdx(current_atom_idx)

            if atom.GetAtomicNum() != 6:
                continue  # Only consider carbon atoms

            neighbors = [nb for nb in atom.GetNeighbors() if nb.GetIdx() != previous_atom_idx]

            # Exclude ester oxygen and any heteroatoms (to ensure unbranched carbon chain)
            for nb in neighbors:
                nb_idx = nb.GetIdx()
                nb_atom_num = nb.GetAtomicNum()

                if nb_idx == ester_o_idx:
                    continue  # Skip ester oxygen
                elif nb_atom_num == 6:
                    stack.append((nb_idx, current_atom_idx))
                else:
                    return False, "Acyl chain is substituted or contains heteroatoms"

            # Limit the chain to prevent infinite loops in cyclic structures
            if len(acyl_chain_atoms) > 8:
                return False, "Acyl chain longer than 8 carbons"

        # Check if the acyl chain has exactly 8 carbon atoms (including carbonyl carbon)
        if len(acyl_chain_atoms) == 8:
            return True, "Contains octanoyl ester group"
        else:
            continue  # Check the next ester group

    return False, "No octanoyl ester groups found"