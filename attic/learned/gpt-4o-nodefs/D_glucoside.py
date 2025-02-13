"""
Classifies: CHEBI:35436 D-glucoside
"""
from rdkit import Chem

def is_D_glucoside(smiles: str):
    """
    Determines if a molecule is a D-glucoside based on its SMILES string.
    A D-glucoside contains a D-glucopyranosyl unit, which can be either 
    alpha or beta, identified by specific stereochemistry.

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

    # Define a SMARTS pattern for both alpha and beta D-glucopyranosides
    # Matching the 6-membered pyranose ring with appropriate stereochemistry
    # Capturing hydroxyl groups and flexible attachment site
    glucoside_pattern = Chem.MolFromSmarts("O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@H](O)[C@H]1O")

    if mol.HasSubstructMatch(glucoside_pattern):
        return True, "Contains D-glucopyranosyl unit indicative of D-glucoside"
    
    return False, "D-glucopyranosyl unit not found"

# Example test case (not run as part of the function definition)
example_smiles = "CC(=O)N[C@@H](CO)[C@@H](OC1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1NC(C)=O)[C@@H](O)[C@H](O)CO"
print(is_D_glucoside(example_smiles))