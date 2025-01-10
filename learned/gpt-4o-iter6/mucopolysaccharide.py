"""
Classifies: CHEBI:37395 mucopolysaccharide
"""
from rdkit import Chem

def is_mucopolysaccharide(smiles: str):
    """ 
    Determines if a molecule is a mucopolysaccharide based on its SMILES string.
    A mucopolysaccharide is a polysaccharide composed of alternating units from uronic acids
    and glycosamines, commonly partially esterified with sulfuric acid.

    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule is a mucopolysaccharide, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (False, "Invalid SMILES string")
    
    # Flexible pattern for uronic acids (sugar with COOH group)
    uronic_acid_pattern = Chem.MolFromSmarts("C1CO[C@@H](CO)[C@H](O)[C@H]1C(=O)[O,H]")
    
    # Improved pattern for glycosamines (sugar with NH2 group)
    glycosamine_pattern = Chem.MolFromSmarts("NC1[C@H](O)[C@H](O)[C@@H](O)[C@H](O1)[C,,O]")

    # Broad pattern for sulfate groups in polysaccharides
    sulfate_pattern = Chem.MolFromSmarts("O[S](=O)(=O)[O]")

    # Assess if molecule contains both uronic acids and glycosamine units
    has_uronic_acid = mol.HasSubstructMatch(uronic_acid_pattern)
    has_glycosamine = mol.HasSubstructMatch(glycosamine_pattern)
    has_sulfate = mol.HasSubstructMatch(sulfate_pattern)

    if has_uronic_acid and has_glycosamine:
        if has_sulfate:
            return (True, "Contains uronic acids and glycosamines with sulfate esterification")
        return (True, "Contains uronic acids and glycosamines but no sulfate esterification")
    else:
        return (False, "No sufficient pattern of uronic acids and glycosamines identified")