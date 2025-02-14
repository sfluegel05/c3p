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
    coa_smarts = "[CX3](=[OX1])S[C@H](O)C(C)(C)COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(=O)(O)O)n1cnc2c(N)ncnc12"
    coa_pattern = Chem.MolFromSmarts(coa_smarts)
    
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA moiety with thioester not found"

    # Define the 3-oxo group SMARTS, including the thioester sulfur.
    oxo_smarts = "[CX3](=[OX1])[CH2][CX3](=[OX1])S"
    oxo_pattern = Chem.MolFromSmarts(oxo_smarts)
    oxo_matches = mol.GetSubstructMatches(oxo_pattern)

    if len(oxo_matches) != 1:
       return False, "Incorrect number of 3-oxo groups connected to thioester found, should be 1"


    # Check if a long carbon chain is attached to oxo_c1_idx
    oxo_c1_idx = oxo_matches[0][0]
    found_chain = False
    for neighbor in mol.GetAtomWithIdx(oxo_c1_idx).GetNeighbors():
        if neighbor.GetSymbol() == 'C':
          found_chain = True
          break
    if not found_chain:
        return False, "No carbon chain attached to 3-oxo group"

    return True, "Contains CoA moiety, 3-oxo group, and a fatty acyl chain"