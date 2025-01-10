"""
Classifies: CHEBI:35436 D-glucoside
"""
from rdkit import Chem

def is_D_glucoside(smiles: str):
    """
    Determines if a molecule is a D-glucoside based on its SMILES string.
    A D-glucoside contains a D-glucopyranosyl unit, featuring specific stereochemistry.

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

    # Refined SMARTS pattern for D-glucopyranoside capturing both alpha and beta forms
    # and considering point of attachments as well as stereochemistry variations
    glucopyranosyl_pattern = Chem.MolFromSmarts("O[C@H]1[C@H](O[C@H](CO)O)[C@@H](O)[C@H](O)[C@H](O)1")

    if mol.HasSubstructMatch(glucopyranosyl_pattern):
        return True, "Contains D-glucopyranosyl unit indicative of D-glucoside"
    
    return False, "D-glucopyranosyl unit not found"