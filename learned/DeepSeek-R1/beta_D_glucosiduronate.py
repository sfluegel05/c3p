"""
Classifies: CHEBI:83411 beta-D-glucosiduronate
"""
"""
Classifies: CHEBI: beta-D-glucosiduronate
"""
from rdkit import Chem

def is_beta_D_glucosiduronate(smiles: str):
    """
    Determines if a molecule is a beta-D-glucosiduronate based on its SMILES string.
    A beta-D-glucosiduronate has a beta-D-glucuronic acid moiety (deprotonated carboxylate) 
    linked via an O-glycosidic bond.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a beta-D-glucosiduronate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Define the beta-D-glucuronide substructure pattern with correct stereochemistry
    glucuronide_pattern = Chem.MolFromSmarts(
        "[O-]C(=O)[C@H]1O[C@H]([C@H]([C@@H]([C@H]1O)O)O)O"
    )
    if glucuronide_pattern is None:
        return False, "Failed to parse glucuronide pattern"
    
    # Check for presence of the glucuronide group
    if mol.HasSubstructMatch(glucuronide_pattern):
        return True, "Contains beta-D-glucuronide group with carboxylate"
    else:
        return False, "No beta-D-glucuronide substructure found"