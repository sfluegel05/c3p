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
        return False, "Invalid SMILES string"

    # Define generalized pattern for glycosidic linkages (can be expanded as needed)
    glycosidic_pattern = Chem.MolFromSmarts("O[C;R][C;R]O")
    
    # Define pattern for amino sugars including diverse aminomonosaccharide context
    amino_sugar_pattern = Chem.MolFromSmarts("[CX4][NX3;H2,H1]")
    
    # Check for glycosidic-like linkages (proxy for polysaccharide structure)
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_pattern)
    
    # Check for presence of amino sugar groups
    amino_sugar_matches = mol.GetSubstructMatches(amino_sugar_pattern)
    
    # Determine classification based on matches
    if len(glycosidic_matches) > 0 and len(amino_sugar_matches) > 0:
        # Use a more representative metric for proportion of amino sugars
        # Relating amino sugar count to glycosidic bonds as a proxy measure
        if len(amino_sugar_matches) / max(1, len(glycosidic_matches)) > 0.2:  # Arbitrary threshold
            return True, "Contains glycosidic-like and substantial proportion of amino sugar-like substructures"
        else:
            return False, "Insufficient proportion of amino sugar-like substructures"
    elif len(glycosidic_matches) == 0:
        return False, "Does not have glycosidic-like substructures"
    elif len(amino_sugar_matches) == 0:
        return False, "Does not contain amino sugar-like substructures"

    return None, "Could not be definitively classified as a glycosaminoglycan"