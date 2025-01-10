"""
Classifies: CHEBI:192499 anthoxanthin
"""
from rdkit import Chem

def is_anthoxanthin(smiles: str):
    """
    Determines if a molecule is an anthoxanthin based on its SMILES string.
    Anthoxanthins are flavonoids characterized by a benzopyranone core structure,
    often with multiple hydroxyl/methoxy groups and potential glycosylation.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an anthoxanthin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Broadened anthoxanthin benzopyranone backbone pattern
    anthoxanthin_pattern = Chem.MolFromSmarts("O=C1C=COC2=CC=CCC12") # Earlier pattern was too restrictive
    
    # Check for backbone structure match
    if not mol.HasSubstructMatch(anthoxanthin_pattern):
        return False, "Does not contain typical anthoxanthin benzopyranone backbone structure"
         
    # Look for hydroxyl or methoxy groups
    hydroxyl_methoxy_pattern = Chem.MolFromSmarts("[OX2H] | [OX1CH3]")  # Focus on OH and OCH3 groups 
    hydroxyl_methoxy_matches = mol.GetSubstructMatches(hydroxyl_methoxy_pattern)
    if len(hydroxyl_methoxy_matches) < 2:  # Typically multiple groups are expected
        return False, "Insufficient hydroxyl or methoxy groups present"
    
    # Check for possible complex glycosylation
    glycoside_pattern = Chem.MolFromSmarts("O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@H]1O")  # Possible glycoside linkages
    if mol.HasSubstructMatch(glycoside_pattern):
        return True, "Anthoxanthin with glycosidic linkage"
    
    # If all checks pass, classify as anthoxanthin
    return True, "Contains anthoxanthin benzopyranone backbone and relevant functional groups"