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

    # Check if the sulfur of the thioester is directly connected to the CoA
    coa_match = mol.GetSubstructMatch(coa_pattern)
    is_connected = False
    
    for match in thioester_matches:
        thioester_sulfur_atom = mol.GetAtomWithIdx(match[1]) # Get the S atom of the thioester
        for neighbor in thioester_sulfur_atom.GetNeighbors(): # Iterate over the S neighbor atoms
           if neighbor.GetIdx() in coa_match: # Check if it belongs to the CoA
               is_connected = True
               break # break the loop if we have found a connected thioester
        if is_connected:
           break # break outer loop if we found one
            
    if not is_connected:
        return False, "Thioester not connected to the CoA moiety"


    return True, "Molecule contains Coenzyme A and a thioester connected to a carboxylic acid."