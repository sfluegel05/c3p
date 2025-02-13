"""
Classifies: CHEBI:18085 glycosaminoglycan
"""
from rdkit import Chem

def is_glycosaminoglycan(smiles: str):
    """
    Determines if a molecule can be classified as a potential glycosaminoglycan 
    based on its SMILES string. Glycosaminoglycans are polysaccharides with a 
    substantial proportion of aminomonosaccharide residues.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a potential glycosaminoglycan, False otherwise
        str: Reason for classification or failure
    """
    
    # Parse SMILES into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Define a simple pattern for polysaccharides 
    polysaccharide_pattern = Chem.MolFromSmarts("O[C@@H]1CO[C@@H](O)[C@H](O)[C@@H]1")
    # Define a simple pattern for amino sugars (e.g., C-N bond)
    amino_sugar_pattern = Chem.MolFromSmarts("[C;!R][NX3;H2,H1]")

    # Check for polysaccharide-like structures using pattern matching
    polysaccharide_match = mol.HasSubstructMatch(polysaccharide_pattern)
    
    # Check for presence of amino groups
    amino_sugar_match = mol.HasSubstructMatch(amino_sugar_pattern)
    
    # Determine classification based on matches
    if polysaccharide_match and amino_sugar_match:
        return True, "Contains polysaccharide-like and amino sugar-like substructures"
    elif not polysaccharide_match:
        return False, "Does not have polysaccharide-like substructures"
    elif not amino_sugar_match:
        return False, "Does not contain amino sugar-like substructures"

    return None, "Could not be definitively classified as a glycosaminoglycan"