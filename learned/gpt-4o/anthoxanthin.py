"""
Classifies: CHEBI:192499 anthoxanthin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_anthoxanthin(smiles: str):
    """
    Determines if a molecule is likely an anthoxanthin based on its SMILES string.
    Anthoxanthins are flavonoid pigments found in plants characterized by a flavonoid core
    and multiple hydroxyl, methoxy, or glycosyl groups.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is likely an anthoxanthin, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Use broader patterns to identify the flavonoid core structure
    flavonoid_patterns = [
        Chem.MolFromSmarts('c1cc(ccc1O)Oc2cc(cc(=O)c2)c3ccccc3'),  # Common flavone pattern
        Chem.MolFromSmarts('c1ccc2c(c1)oc(=O)c3cc(O)cc(O)c23')     # Chromone/flavone pattern
    ]
    
    if not any(mol.HasSubstructMatch(pattern) for pattern in flavonoid_patterns):
        return False, "No flavonoid core structure found"
    
    # Search for the presence of at least two hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts('[OX2H]')
    num_hydroxyls = len(mol.GetSubstructMatches(hydroxyl_pattern))
    if num_hydroxyls < 2:
        return False, f"Insufficient number of hydroxy groups, found {num_hydroxyls}"
    
    # Search for methoxy groups
    methoxy_pattern = Chem.MolFromSmarts('COC')
    num_methoxy = len(mol.GetSubstructMatches(methoxy_pattern))
    
    # Search for glycosidic linkages indicative of glycosylated structures
    glycoside_pattern = Chem.MolFromSmarts('[C@H]1([O-])[C@@H](O1)')
    num_glycosides = len(mol.GetSubstructMatches(glycoside_pattern))
    
    # Use flexible logic to assess functionality common to anthoxanthins
    if num_glycosides > 0 or num_methoxy > 0:
        return True, "Contains flavonoid core with hydroxy, methoxy or glycosyl groups, indicating an anthoxanthin"

    # Consider the criteria for molecules with multiple functional groups
    if num_hydroxyls >= 2 and (num_methoxy > 0 or num_glycosides > 0):
        return True, "Contains flavonoid core with several functional groups typical of anthoxanthins"

    return False, "Lacks sufficient functional group diversity for classification as anthoxanthin"