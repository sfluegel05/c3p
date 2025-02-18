"""
Classifies: CHEBI:47908 alkanethiol
"""
"""
Classifies: CHEBI:46643 alkanethiol
"""
from rdkit import Chem

def is_alkanethiol(smiles: str):
    """
    Determines if a molecule is an alkanethiol based on its SMILES string.
    An alkanethiol has at least one sulfanyl group (-SH) attached to an alkyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkanethiol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Find all SH groups
    sh_pattern = Chem.MolFromSmarts('[SH]')
    matches = mol.GetSubstructMatches(sh_pattern)
    
    for s_idx, _ in matches:
        s_atom = mol.GetAtomWithIdx(s_idx)
        # Check if S is bonded to exactly one H and one C
        if s_atom.GetDegree() != 2:
            continue
        h_count = sum(1 for neighbor in s_atom.GetNeighbors() if neighbor.GetAtomicNum() == 1)
        if h_count != 1:
            continue
        
        # Get the attached carbon
        carbon = None
        for neighbor in s_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:
                carbon = neighbor
                break
        if not carbon:
            continue
        
        # Check carbon is part of an alkyl group (all adjacent atoms are C or H except for the S)
        valid = True
        for neighbor in carbon.GetNeighbors():
            if neighbor.GetIdx() == s_idx:
                continue
            if neighbor.GetAtomicNum() not in [6, 1]:
                valid = False
                break
        if valid:
            return True, "Contains a sulfanyl group (-SH) attached to an alkyl group"
    
    return False, "No sulfanyl group attached to an alkyl group"