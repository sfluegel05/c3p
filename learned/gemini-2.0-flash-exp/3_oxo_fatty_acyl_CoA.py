"""
Classifies: CHEBI:15489 3-oxo-fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
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

    # Define the CoA moiety SMARTS (simplified, focusing on key linkage)
    coa_smarts = "[C@H](O)C(C)(C)COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(=O)(O)O)n1cnc2c(N)ncnc12"
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
    thioester_smarts = "C(=O)S"
    thioester_pattern = Chem.MolFromSmarts(thioester_smarts)
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)

    if len(thioester_matches) != 1:
         return False, "Incorrect number of thioester groups, should be 1"

    # Check the connection between 3-oxo, thioester, and CoA
    oxo_c1_idx = oxo_matches[0][0]
    oxo_c2_idx = oxo_matches[0][2]
    thioester_c_idx = thioester_matches[0][0]
    thioester_s_idx = thioester_matches[0][1]
    
    if thioester_c_idx != oxo_c2_idx:
      return False, "Thioester not connected to the 3-oxo group"

    # Check if a long carbon chain is attached to oxo_c1_idx
    
    found_chain = False
    for neighbor in mol.GetAtomWithIdx(oxo_c1_idx).GetNeighbors():
        if neighbor.GetSymbol() == 'C':
          found_chain = True
          break
    if not found_chain:
        return False, "No carbon chain attached to 3-oxo group"
    

    # Verify that the sulfur of the thioester is connected to the pantethiene of the CoA
    found_coa_connection = False
    for coa_match in mol.GetSubstructMatches(coa_pattern):
        for atom_idx in coa_match:
          if mol.GetAtomWithIdx(atom_idx).GetSymbol() == 'S':
                # Check for direct connection to the thioester S
                for neighbor_idx in mol.GetAtomWithIdx(atom_idx).GetNeighbors():
                    if neighbor_idx == thioester_s_idx:
                        found_coa_connection = True
                        break
                if found_coa_connection:
                    break #no need to check the other atoms of the CoA moiety

        if found_coa_connection:
            break #no need to check the other matches to CoA


    if not found_coa_connection:
      return False, "Thioester and CoA not connected"
    
    # Count rotatable bonds to verify long chains, this time is more relaxed
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Chains too short to be fatty acids"


    # Check molecular weight - acyl-CoA typically >500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
      return False, "Molecular weight too low for acyl-CoA"

    return True, "Contains CoA moiety, 3-oxo group, and a fatty acyl chain"