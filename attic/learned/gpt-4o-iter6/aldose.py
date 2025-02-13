"""
Classifies: CHEBI:15693 aldose
"""
from rdkit import Chem

def is_aldose(smiles: str):
    """
    Determines if a molecule is an aldose based on its SMILES string.
    An aldose is a polyhydroxy aldehyde or its intramolecular hemiacetal.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an aldose, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Linear aldehyde pattern - terminal C(=O)H
    aldehyde_pattern = Chem.MolFromSmarts("[#6][CX3H](=O)[#6]")
    
    # Cyclic hemiacetal pattern - allowing 5 or 6-membered rings with oxygen
    furanose_pattern = Chem.MolFromSmarts("O1CCCC1")
    pyranose_pattern = Chem.MolFromSmarts("O1CCCCC1")

    # Check for the presence of an aldehyde group or furanose/pyranose ring structure
    if not (mol.HasSubstructMatch(aldehyde_pattern) or 
            mol.HasSubstructMatch(furanose_pattern) or 
            mol.HasSubstructMatch(pyranose_pattern)): 
        return False, "No aldehyde group or cyclic form detected"

    # Check for sufficient hydroxyl groups (polyhydroxy structure)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 2:
        return False, "Insufficient number of hydroxyl groups for polyhydroxy structure"

    return True, "Contains structural features consistent with aldose - polyhydroxy aldehyde or cyclic form"