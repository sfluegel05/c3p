"""
Classifies: CHEBI:50753 isoflavonoid
"""
from rdkit import Chem

def is_isoflavonoid(smiles: str):
    """
    Determines if a molecule is an isoflavonoid based on its SMILES string.
    An isoflavonoid is defined as a 1-benzopyran with an aryl substituent at position 3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an isoflavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generalized 1-benzopyran structure, looking for a bicyclic system with oxygen in one ring
    benzopyran_pattern = Chem.MolFromSmarts("c1cc2occc2c1")
    if benzopyran_pattern is None:
        return False, "Error in SMARTS pattern definition"

    if not mol.HasSubstructMatch(benzopyran_pattern):
        return False, "No 1-benzopyran backbone found"
        
    # Generic pattern for an aryl group (a phenyl ring) attached
    aryl_group_pattern = Chem.MolFromSmarts("c1ccccc1")
    if aryl_group_pattern is None:
        return False, "Error in SMARTS pattern definition for aryl group"

    if not mol.HasSubstructMatch(aryl_group_pattern):
        return False, "No aryl substituent at position 3 found"

    # Additional check to ensure aryl group is at the correct position (3)
    # Using a smarter substructure search strategy
    for match in mol.GetSubstructMatches(benzopyran_pattern):
        # Identify the part of the molecule that matches the benzopyran
        # and check adjacency of the aryl group
        # Assuming the second atom in the match is the oxygen atom in the benzopyran
        if len(match) >= 3 and mol.GetAtomWithIdx(match[3]).HasSubstructMatch(aryl_group_pattern):
            return True, "Contains 1-benzopyran with an aryl substituent at position 3"
    
    return False, "Correct aryl substituent placement not confirmed"

# Testing a sample isoflavonoid SMILES
result, reason = is_isoflavonoid("COc1ccc(cc1)-c1coc2cc(O)cc(O)c2c1=O")
print(result, reason)