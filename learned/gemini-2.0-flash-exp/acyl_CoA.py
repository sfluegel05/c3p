"""
Classifies: CHEBI:17984 acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_acyl_CoA(smiles: str):
    """
    Determines if a molecule is an acyl-CoA based on its SMILES string.
    An acyl-CoA is a thioester formed between coenzyme A and a carboxylic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for the CoA moiety.
    coa_pattern = Chem.MolFromSmarts("CCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)([OH])OP(=O)([OH])OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(=O)([OH])O)n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Coenzyme A substructure not found"

    # Define a SMARTS pattern for the thioester group (C(=O)-S)
    thioester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[SX2]")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
         return False, "Thioester group not found"

    # Check if the carbonyl carbon of the thioester is connected to the CoA
    coa_match = mol.GetSubstructMatch(coa_pattern)
    is_connected = False

    for match in thioester_matches:
       carbonyl_carbon_idx = match[0] # Get the C atom of the C=O group
       thioester_sulfur_idx = match[1] # Get the S atom of the C-S group
       carbonyl_carbon_atom = mol.GetAtomWithIdx(carbonyl_carbon_idx)
       sulfur_atom = mol.GetAtomWithIdx(thioester_sulfur_idx)

       if carbonyl_carbon_atom.GetIdx() in coa_match:
          for neighbor in carbonyl_carbon_atom.GetNeighbors():
              if neighbor.GetIdx() == sulfur_atom.GetIdx():
                is_connected = True
                break # if we found the sulfur break
       if is_connected:
          break  # break out of outer loop


    if not is_connected:
        return False, "Thioester not connected to the CoA moiety"


    return True, "Molecule contains Coenzyme A and a thioester connected to a carboxylic acid."