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
    
    # Parse the SMILES string to a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # More flexible aldehyde pattern
    aldehyde_pattern = Chem.MolFromSmarts("[CH]=O")
    
    # Updated furanose pattern considering oxygen linkage
    furanose_pattern = Chem.MolFromSmarts("OC1C[C@H](O)[C@@H]O1")
    pyranose_pattern = Chem.MolFromSmarts("OC1C[C@H](O)[C@@H](O)C1O")

    # Check for the presence of an aldehyde group or modified furanose/pyranose ring structure
    if not (mol.HasSubstructMatch(aldehyde_pattern) or 
            mol.HasSubstructMatch(furanose_pattern) or 
            mol.HasSubstructMatch(pyranose_pattern)): 
        return False, "No aldehyde group or cyclic form detected"

    # Check for sufficient hydroxyl groups (polyhydroxy structure)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 3:
        return False, f"Insufficient hydroxyl groups for polyhydroxy structure, found {len(hydroxyl_matches)}"

    return True, "Contains structural features consistent with aldose - polyhydroxy aldehyde or cyclic form"