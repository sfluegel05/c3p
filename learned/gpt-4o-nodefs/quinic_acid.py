"""
Classifies: CHEBI:26493 quinic acid
"""
from rdkit import Chem

def is_quinic_acid(smiles: str):
    """
    Determines if a molecule is a quinic acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a quinic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for cyclohexane core with multiple hydroxyls and one carboxylic acid
    # Quinic acid base structure: OC1CC(O)(C(CCC1O)C(O)=O)O
    quinic_acid_base_pattern = Chem.MolFromSmarts("OC1CCC(O)(C(O)=O)CC1O")
    
    if not mol.HasSubstructMatch(quinic_acid_base_pattern):
        return False, "Mismatch in quinic acid core structure (cyclohexane with hydroxyls and carboxylic acid)"

    # If matches, check for any additional features or derivatives
    # Note: This part is flexible as quinic acid derivatives can vary
    # For simple quinic acids
    base_match = mol.HasSubstructMatch(quinic_acid_base_pattern)
    if base_match:
        return True, "Matches quinic acid core structure"

    return False, "Does not match essential quinic acid features"