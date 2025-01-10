"""
Classifies: CHEBI:35436 D-glucoside
"""
from rdkit import Chem

def is_D_glucoside(smiles: str):
    """
    Determines if a molecule is a D-glucoside based on its SMILES string.
    A D-glucoside typically contains a D-glucopyranosyl unit that may appear
    as a specific substructure in a larger molecule, often with specific stereochemistry.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a D-glucoside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a more specific and flexible SMARTS pattern for D-glucopyranosyl rings
    # Capturing both alpha and beta anomers with a common D-glucose pattern
    # The below SMARTS attempts to match the glucose unit with variable substituents
    glucoside_pattern = Chem.MolFromSmarts("[O;R][C@H]1[C@@H]([C@H]([C@@H](O)[C@H](O)[C@H]1O)O)(*[C@H](C)*)")
    
    if mol.HasSubstructMatch(glucoside_pattern):
        return True, "Contains D-glucopyranosyl unit indicative of D-glucoside"
    
    return False, "D-glucopyranosyl unit not found"

# Test cases should include manual validation to ensure diverse structures are correctly identified