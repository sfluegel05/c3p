"""
Classifies: CHEBI:18085 glycosaminoglycan
"""
from rdkit import Chem

def is_glycosaminoglycan(smiles: str):
    """
    Determines if a molecule is a glycosaminoglycan based on its SMILES string.
    A glycosaminoglycan is characterized by polysaccharides containing substantial proportions of aminomonosaccharide residues.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a glycosaminoglycan, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify sugar-like patterns
    sugar_pattern = Chem.MolFromSmarts("C(O)C(O)C")
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No polysaccharide pattern detected"

    # Identify amino sugar groups
    amino_sugar_pattern = Chem.MolFromSmarts("[NX3;H2]")  # Simple NH2 pattern
    amino_sugar_matches = mol.GetSubstructMatches(amino_sugar_pattern)
    
    if len(amino_sugar_matches) < 1:
        return False, "No aminomonosaccharide residues found"

    # Assuming multiple sugar units and amino groups qualify the molecule
    # However, we might need to increase the threshold based on requirements
    if len(amino_sugar_matches) >= 2:  # Threshold can be adjusted
        return True, "Contains polysaccharide chain with aminomonosaccharide residues"
    
    return False, "Does not meet the threshold for aminomonosaccharide residues"

# Example usage
# status, reason = is_glycosaminoglycan("your_SMILES_here")