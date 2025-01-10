"""
Classifies: CHEBI:47622 acetate ester
"""
"""
Classifies: CHEBI:33282 acetate ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_acetate_ester(smiles: str):
    """
    Determines if a molecule is an acetate ester based on its SMILES string.
    An acetate ester contains the acetate group (-OC(=O)CH3) attached to a carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an acetate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for acetate ester pattern: -OC(=O)CH3
    # [OX2] oxygen with 2 connections
    # [CX3] carbon with 3 connections (sp2)
    # [=OX1] oxygen double bond
    # [CH3] methyl group
    acetate_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])[CH3]")
    matches = mol.GetSubstructMatches(acetate_pattern)
    
    if not matches:
        return False, "No acetate ester group found"
    
    # For each match, verify it's a proper acetate ester
    valid_acetates = 0
    for match in matches:
        o_atom = mol.GetAtomWithIdx(match[0])
        c_atom = mol.GetAtomWithIdx(match[1])
        
        # Check if oxygen is connected to carbon (not hydrogen or other atoms)
        for neighbor in o_atom.GetNeighbors():
            if neighbor.GetIdx() != match[1]:  # not the carbonyl carbon
                if neighbor.GetAtomicNum() == 6:  # is carbon
                    valid_acetates += 1
                    break
    
    if valid_acetates == 0:
        return False, "Acetate group present but not properly connected"
        
    # Additional check to exclude acetic acid itself
    if len(mol.GetAtoms()) <= 4:  # C2H4O2 (acetic acid) has 4 atoms
        return False, "Molecule is too small - might be acetic acid"
        
    # Check for carboxylic acid pattern to ensure we're not counting it
    acid_pattern = Chem.MolFromSmarts("[OH][CX3](=[OX1])[#6]")
    if mol.HasSubstructMatch(acid_pattern):
        if len(mol.GetSubstructMatches(acid_pattern)) == len(matches):
            return False, "Contains carboxylic acid group(s) but no acetate ester"
            
    return True, f"Contains {valid_acetates} acetate ester group(s)"