"""
Classifies: CHEBI:83411 beta-D-glucosiduronate
"""
"""
Classifies: beta-D-glucosiduronate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_beta_D_glucosiduronate(smiles: str):
    """
    Determines if a molecule contains a beta-D-glucosiduronate group.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule contains beta-D-glucosiduronate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for beta-D-glucuronic acid
    # Pattern matches:
    # - O-glycosidic bond
    # - Beta configuration at C1
    # - Correct stereochemistry at C2,C3,C4
    # - Carboxylate group at C5
    glucuronate_pattern = Chem.MolFromSmarts(
        "[OX2][C@@H]1[C@@H]([C@H]([C@@H]([C@H](O1)C(=O)[O-])O)O)O"
    )
    
    if not mol.HasSubstructMatch(glucuronate_pattern):
        return False, "No beta-D-glucuronate group found with correct stereochemistry"
    
    # Count matches to ensure we have at least one
    matches = mol.GetSubstructMatches(glucuronate_pattern)
    if len(matches) < 1:
        return False, "No beta-D-glucuronate group found"
        
    # Verify presence of carboxylate
    carboxylate_pattern = Chem.MolFromSmarts("C(=O)[O-]")
    if not mol.HasSubstructMatch(carboxylate_pattern):
        return False, "No carboxylate group found"
    
    # Check for required hydroxyls
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    hydroxyls = len(mol.GetSubstructMatches(hydroxyl_pattern))
    if hydroxyls < 3:
        return False, "Insufficient hydroxyl groups"
        
    # Verify glycosidic bond
    glycosidic_pattern = Chem.MolFromSmarts("[OX2][CH]1O[CH][CH][CH][CH]1")
    if not mol.HasSubstructMatch(glycosidic_pattern):
        return False, "No glycosidic bond found"

    return True, "Contains beta-D-glucuronate group with correct stereochemistry and carboxylate"