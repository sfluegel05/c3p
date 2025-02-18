"""
Classifies: CHEBI:22798 beta-D-glucoside
"""
"""
Classifies: CHEBI:XXXXX beta-D-glucoside
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

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
    # Matches beta-D-glucopyranose core with correct stereochemistry
    # Allows substituents on hydroxyl groups (e.g., acetylations)
    beta_D_glucoside_pattern = Chem.MolFromSmarts(
        "[OX2;!H0]-[C@@H]1[C@@H]([OX2])[C@H]([OX2])[C@@H]([OX2])[C@H](C[OX2])O1"
    )
    
    if beta_D_glucoside_pattern is None:
        return False, "Invalid SMARTS pattern"
    
    # Check for the presence of the beta-D-glucoside pattern
    if mol.HasSubstructMatch(beta_D_glucoside_pattern):
        return True, "Beta-D-glucoside detected: anomeric carbon in beta configuration with D-glucose stereochemistry"
    
    # Additional check for possible modifications (e.g., C6 as CH2OR)
    # This pattern allows any substituent on the C6 oxygen
    # Ensure the glucose stereochemistry is maintained
    # Alternative pattern to capture different substitution scenarios
    alt_pattern = Chem.MolFromSmarts(
        "[OX2;!H0]-[C@@H]1[C@@H]([!#1])[C@H]([!#1])[C@@H]([!#1])[C@H](C[OX2])O1"
    )
    if alt_pattern and mol.HasSubstructMatch(alt_pattern):
        return True, "Beta-D-glucoside detected: anomeric carbon in beta configuration with D-glucose stereochemistry"
    
    return False, "No beta-D-glucoside substructure found"