"""
Classifies: CHEBI:37395 mucopolysaccharide
"""
from rdkit import Chem

def is_mucopolysaccharide(smiles: str):
    """ 
    Determines if a molecule is a mucopolysaccharide based on SMILES string.
    A mucopolysaccharide is a polysaccharide composed of alternating units from uronic acids
    and glycosamines, commonly partially esterified with sulfuric acid.

    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a mucopolysaccharide, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Patterns for detecting uronic acids
    # Carbohydrate ring with a carboxylic acid (COOH) group
    uronic_acid_pattern = Chem.MolFromSmarts("C1([C@@H]([C@H]([C@H](O1)CO)O)O)C(=O)[O]")
    
    # Patterns for detecting glycosamine
    # Nitrogen attached in sugar-like ring structure
    glycosamine_pattern = Chem.MolFromSmarts("[C@@H]1([C@H]([C@H]([C@H](C(O1)O)NC(=O))O)O)")

    # Sulfate group esterified in a polysaccharide chain
    sulfate_pattern = Chem.MolFromSmarts("OS(=O)(=O)[O,C]")

    # Assess if molecule contains both uronic acids and glycosamines
    has_uronic_acid = mol.HasSubstructMatch(uronic_acid_pattern)
    has_glycosamine = mol.HasSubstructMatch(glycosamine_pattern)
    has_sulfate = mol.HasSubstructMatch(sulfate_pattern)

    if has_uronic_acid and has_glycosamine:
        if has_sulfate:
            return True, "Contains uronic acids and glycosamines with sulfate esterification"
        return True, "Contains uronic acids and glycosamines but no sulfate esterification"
    else:
        return False, "No sufficient pattern of uronic acids and glycosamines identified"