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
        bool: True if the molecule is a phenyl acetate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define ester SMARTS pattern
    ester_smarts = '[CX3](=O)[OX2H0]'
    ester_pattern = Chem.MolFromSmarts(ester_smarts)
    if ester_pattern is None:
        return False, "Invalid ester SMARTS pattern"

    # Search for ester groups in the molecule
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester groups found"

    # Iterate over ester groups
    for match in ester_matches:
        carbonyl_c_idx = match[0]
        ester_o_idx = match[1]

        carbonyl_c = mol.GetAtomWithIdx(carbonyl_c_idx)
        ester_o = mol.GetAtomWithIdx(ester_o_idx)

        # Check if acyl group is acetyl (carbonyl carbon connected to methyl group)
        is_acetyl = False
        for neighbor in carbonyl_c.GetNeighbors():
            if neighbor.GetIdx() == ester_o_idx:
                continue
            if neighbor.GetAtomicNum() == 6 and neighbor.GetDegree() == 1 and neighbor.GetTotalNumHs() == 3:
                is_acetyl = True
                break
        if not is_acetyl:
            continue  # Not an acetyl group

        # Check if ester oxygen is connected to an aromatic carbon
        for neighbor in ester_o.GetNeighbors():
            if neighbor.GetIdx() == carbonyl_c_idx:
                continue
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIsAromatic():
                return True, "Contains phenyl acetate group"

    return False, "No phenyl acetate group found"