"""
Classifies: CHEBI:33447 phospho sugar
"""
"""
Classifies: CHEBI:phospho sugar
"""
from rdkit import Chem
from rdkit.Chem import MolFromSmiles, MolFromSmarts

def is_phospho_sugar(smiles: str):
    """
    Determines if a molecule is a phospho sugar (monosaccharide with phosphate ester).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if phospho sugar, False otherwise
        str: Reason for classification
    """
    mol = MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES"
    
    # Check for monosaccharide: cyclic structure with multiple hydroxyls
    # Using a generic sugar pattern: 5 or 6-membered ring with multiple O's and OH groups
    sugar_pattern = MolFromSmarts("[O;R]@;=,#[C;R]@;=,#[C;R]@;=,#[C;R]@;=,#[C;R]@;=,#[O;R]")
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No sugar ring detected"
    
    # Check for phosphate ester: O-P(=O)(O)-O connected to carbon
    phosphate_ester = MolFromSmarts("[C!$(*=O)][OX2][P](=[OX1])([OX2-])[OX2]")
    if not mol.HasSubstructMatch(phosphate_ester):
        return False, "No phosphate ester group found"
    
    # Ensure phosphate is attached to an alcoholic oxygen (not carboxylic acid, etc.)
    # Check that the oxygen connected to P is bonded to a carbon without double bonds
    alcoholic_phosphate = MolFromSmarts("[CX4][OX2][P]")
    if not mol.HasSubstructMatch(alcoholic_phosphate):
        return False, "Phosphate not on alcoholic hydroxyl"
    
    return True, "Monosaccharide with phosphate ester"