"""
Classifies: CHEBI:35436 D-glucoside
"""
from rdkit import Chem

def is_D_glucoside(smiles: str):
    """
    Determines if a molecule is a D-glucoside based on its SMILES string.
    A D-glucoside typically contains a D-glucopyranosyl unit that may appear
    as a specific substructure in a larger molecule.
    
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
    
    # Define the D-glucopyranoside substructure pattern
    # This pattern needs to target the specific glucose ring with correct stereochemistry
    # including beta anomeric configuration for D-glucosides
    glucoside_pattern = Chem.MolFromSmarts("O[C@H]1[C@@H]([C@H](O)[C@H](O)[C@H](CO)O1)C")
    
    if mol.HasSubstructMatch(glucoside_pattern):
        return True, "Contains D-glucopyranosyl unit indicative of D-glucoside"
    
    return False, "D-glucopyranosyl unit not found"

# Note: This output assumes SMILES input follows conventional representation
# for D-glucopyranosyl units. Adjustments to SMARTS may be required.