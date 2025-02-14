"""
Classifies: CHEBI:15489 3-oxo-fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3_oxo_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-oxo-fatty acyl-CoA based on its SMILES string.
    A 3-oxo-fatty acyl-CoA has a CoA moiety, a 3-oxo group and a fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the CoA moiety SMARTS, including the thioester sulfur and carbonyl
    coa_smarts = "[C@H](O)C(C)(C)COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(=O)(O)O)n1cnc2c(N)ncnc12C(=O)S"
    coa_pattern = Chem.MolFromSmarts(coa_smarts)
    
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA moiety not found"

    # Define the 3-oxo group SMARTS, and also the carbonyl of the thioester.
    oxo_smarts = "[CX3](=[OX1])[CH2][CX3](=[OX1])"
    oxo_pattern = Chem.MolFromSmarts(oxo_smarts)
    oxo_matches = mol.GetSubstructMatches(oxo_pattern)

    if len(oxo_matches) != 1:
       return False, "Incorrect number of 3-oxo groups found, should be 1"

    # Check for thioester group connected to CoA
    thioester_smarts = "[CX3](=[OX1])S"
    thioester_pattern = Chem.MolFromSmarts(thioester_smarts)
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)

    if len(thioester_matches) != 1:
         return False, "Incorrect number of thioester groups, should be 1"

    # Check the connection between 3-oxo, thioester, and CoA
    oxo_c1_idx = oxo_matches[0][0]
    oxo_c2_idx = oxo_matches[0][1] #methylene carbon
    thioester_c_idx = thioester_matches[0][0]

    found_correct_connection = False
    for neighbor in mol.GetAtomWithIdx(oxo_c2_idx).GetNeighbors():
        if neighbor.GetIdx() == thioester_c_idx:
            found_correct_connection = True
            break
    if not found_correct_connection:
        return False, "Thioester not connected to the 3-oxo group"

    # Check if a long carbon chain is attached to oxo_c1_idx
    
    found_chain = False
    for neighbor in mol.GetAtomWithIdx(oxo_c1_idx).GetNeighbors():
        if neighbor.GetSymbol() == 'C':
          found_chain = True
          break
    if not found_chain:
        return False, "No carbon chain attached to 3-oxo group"

    return True, "Contains CoA moiety, 3-oxo group, and a fatty acyl chain"