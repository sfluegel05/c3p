"""
Classifies: CHEBI:35759 1-monoglyceride
"""
"""
Classifies: CHEBI:17854 1-monoglyceride
"""

from rdkit import Chem

def is_1_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 1-monoglyceride based on its SMILES string.
    A 1-monoglyceride is a monoglyceride in which the acyl substituent is located at position 1.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-monoglyceride, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the ester functional group SMARTS pattern
    ester_smarts = '[#6][C](=O)[O][CH2]'
    ester_pattern = Chem.MolFromSmarts(ester_smarts)
    ester_matches = mol.GetSubstructMatches(ester_pattern)

    if not ester_matches:
        return False, "No ester group connected to a primary carbon found"

    found_glycerol = False

    # Iterate over ester matches to find the glycerol backbone
    for match in ester_matches:
        acyl_chain_c = match[0]
        carbonyl_c = match[1]
        ester_o = match[2]
        glycerol_c1 = match[3]  # Primary carbon (position 1)

        glycerol_c1_atom = mol.GetAtomWithIdx(glycerol_c1)
        # Check if glycerol_c1 is connected to exactly one other carbon (position 2)
        neighbors_c1 = [nbr.GetIdx() for nbr in glycerol_c1_atom.GetNeighbors() if nbr.GetIdx() != ester_o]
        if len(neighbors_c1) != 1:
            continue
        c2_idx = neighbors_c1[0]
        c2_atom = mol.GetAtomWithIdx(c2_idx)
        if c2_atom.GetAtomicNum() != 6:
            continue

        # Check neighbors of position 2 carbon (should be connected to position 1, position 3, and a hydroxyl group)
        neighbors_c2 = [nbr.GetIdx() for nbr in c2_atom.GetNeighbors() if nbr.GetIdx() != glycerol_c1]
        if len(neighbors_c2) != 2:
            continue

        c3_idx = None
        o2_idx = None
        for idx in neighbors_c2:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 6:
                c3_idx = idx
            elif atom.GetAtomicNum() == 8:
                o2_idx = idx
        if c3_idx is None or o2_idx is None:
            continue

        # Verify that the oxygen at position 2 is a hydroxyl group
        o2_atom = mol.GetAtomWithIdx(o2_idx)
        if o2_atom.GetTotalDegree() != 1 or o2_atom.GetImplicitValence() != 1:
            continue

        # Check neighbors of position 3 carbon (should be connected to position 2 and a hydroxyl group)
        c3_atom = mol.GetAtomWithIdx(c3_idx)
        neighbors_c3 = [nbr.GetIdx() for nbr in c3_atom.GetNeighbors() if nbr.GetIdx() != c2_idx]
        if len(neighbors_c3) != 1:
            continue
        o3_idx = neighbors_c3[0]
        o3_atom = mol.GetAtomWithIdx(o3_idx)
        if o3_atom.GetAtomicNum() != 8:
            continue
        # Verify that the oxygen at position 3 is a hydroxyl group
        if o3_atom.GetTotalDegree() != 1 or o3_atom.GetImplicitValence() != 1:
            continue

        found_glycerol = True
        break

    if found_glycerol:
        return True, "Contains glycerol backbone with acyl group at position 1"
    else:
        return False, "Does not contain glycerol backbone with acyl group at position 1"