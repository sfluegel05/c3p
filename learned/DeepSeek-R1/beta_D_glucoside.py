"""
Classifies: CHEBI:22798 beta-D-glucoside
"""
"""
Classifies: CHEBI:XXXXX beta-D-glucoside
"""
from rdkit import Chem

def is_beta_D_glucoside(smiles: str):
    """
    Determines if a molecule is a beta-D-glucoside based on its SMILES string.
    A beta-D-glucoside has a glucose moiety where the anomeric carbon (C1) has a beta configuration.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a beta-D-glucoside, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Define SMARTS pattern for beta-D-glucoside
    # Pattern matches the anomeric carbon (C1) with beta configuration (C@@) and D-glucose stereochemistry
    beta_D_glucoside_pattern = Chem.MolFromSmarts(
        "[OX2;!H0]-[C@@H]1[C@@H]([!#1])[C@H]([!#1])[C@@H]([!#1])[C@H](CO)O1"
    )
    
    if beta_D_glucoside_pattern is None:
        return False, "Invalid SMARTS pattern"
    
    # Check for the presence of the beta-D-glucoside pattern
    if mol.HasSubstructMatch(beta_D_glucoside_pattern):
        return True, "Beta-D-glucoside detected: anomeric carbon in beta configuration with D-glucose stereochemistry"
    
    return False, "No beta-D-glucoside substructure found"