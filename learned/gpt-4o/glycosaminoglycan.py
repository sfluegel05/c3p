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

    # Define a more comprehensive pattern for polysaccharides (focusing on acetal linkages)
    polysaccharide_pattern = Chem.MolFromSmarts("O[C@@H]1CO[C@@H](O)[C@H](O)[C@@H]1")
    
    # Define a more precise pattern for amino sugars (aminomonosaccharides)
    amino_sugar_pattern = Chem.MolFromSmarts("[CX4H][NX3;H2,H1][CX4H]")  # Slightly more rigorous to amino sugar context

    # Check for polysaccharide-like structures using pattern matching
    polysaccharide_match = mol.HasSubstructMatch(polysaccharide_pattern)
    
    # Check for presence of amino sugar groups
    amino_sugar_matches = mol.GetSubstructMatches(amino_sugar_pattern)
    
    # Determine classification based on matches
    if polysaccharide_match and len(amino_sugar_matches) > 0:
        # Check substantial proportion of amino sugars
        total_osmol = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8])  # Potential over count, just an indicator
        if len(amino_sugar_matches) / total_osmol > 0.1:  # Arbitrary threshold for 'substantial'
            return True, "Contains polysaccharide-like and sufficient proportion of amino sugar-like substructures"
        else:
            return False, "Amino groups present but not in substantial proportion"
    elif not polysaccharide_match:
        return False, "Does not have polysaccharide-like substructures"
    elif len(amino_sugar_matches) == 0:
        return False, "Does not contain amino sugar-like substructures"

    return None, "Could not be definitively classified as a glycosaminoglycan"