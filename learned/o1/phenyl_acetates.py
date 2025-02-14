"""
Classifies: CHEBI:140310 phenyl acetates
"""
"""
Classifies: Phenyl Acetates
"""
from rdkit import Chem

def is_phenyl_acetates(smiles: str):
    """
    Determines if a molecule is a phenyl acetate based on its SMILES string.
    A phenyl acetate is defined as an acetate ester obtained by formal condensation 
    of the carboxy group of acetic acid with the hydroxy group of any phenol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phenyl acetate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify ester functional groups
    ester_smarts = '[#6][C](=O)[O][#6]'  # General ester pattern
    ester_pattern = Chem.MolFromSmarts(ester_smarts)
    ester_matches = mol.GetSubstructMatches(ester_pattern)

    if not ester_matches:
        return False, "No ester groups found"

    # Check each ester group
    for match in ester_matches:
        # Indices of atoms in the ester group
        c1_idx, c2_idx, o_idx, c3_idx = match

        # Get atom objects
        c_acyl = mol.GetAtomWithIdx(c2_idx)  # Carbonyl carbon
        o_ester = mol.GetAtomWithIdx(o_idx)  # Ester oxygen
        c_alkyl = mol.GetAtomWithIdx(c3_idx)  # Alkoxy carbon

        # Check if acyl carbon is part of acetic acid (attached to methyl group)
        acyl_neighbors = [atom for atom in c_acyl.GetNeighbors() if atom.GetIdx() != o_ester.GetIdx()]
        has_methyl_group = False
        for neighbor in acyl_neighbors:
            if neighbor.GetAtomicNum() == 6 and neighbor.GetDegree() == 3:
                # Check if all neighbors are hydrogens (methyl group)
                hydrogen_count = sum(1 for n in neighbor.GetNeighbors() if n.GetAtomicNum() == 1)
                if hydrogen_count == 3:
                    has_methyl_group = True
                    break
        if not has_methyl_group:
            continue  # Not acetic acid

        # Check if alkoxy part is connected to aromatic ring (phenol derivative)
        # Get the atom connected to the ester oxygen besides the carbonyl carbon
        o_neighbors = [atom for atom in o_ester.GetNeighbors() if atom.GetIdx() != c_acyl.GetIdx()]
        if len(o_neighbors) != 1:
            continue  # Unusual ester oxygen connectivity
        alkoxy_atom = o_neighbors[0]

        # Check if alkoxy atom is an aromatic carbon
        if alkoxy_atom.GetAtomicNum() == 6 and alkoxy_atom.GetIsAromatic():
            return True, "Contains phenyl acetate group"

    return False, "Does not match phenyl acetate criteria"