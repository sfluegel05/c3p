"""
Classifies: CHEBI:26493 quinic acid
"""
from rdkit import Chem

def is_quinic_acid(smiles: str):
    """
    Determines if a molecule is quinic acid based on its SMILES string.
    Quinic acid is defined as a cyclitol carboxylic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is quinic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for a cyclohexane ring
    ring_info = mol.GetRingInfo()
    if not any(len(ring) == 6 for ring in ring_info.AtomRings()):
        return False, "No cyclohexane ring found"

    # Look for hydroxy groups on the cyclohexane ring
    hydroxy_pattern = Chem.MolFromSmarts("[CX4;R][OX2H]")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    if len(hydroxy_matches) < 3:  # Quinic acid typically has 3 or more hydroxyls
        return False, f"Not enough hydroxyl groups found, got {len(hydroxy_matches)}"

    # Check for carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O;H1]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # We might want to check additional patterns for esterified quinic acid derivatives

    return True, "Contains features consistent with quinic acid"

# Example SMILES (quinic acid):
smiles_example = "O[C@H]1C[C@@](O)(C[C@H](O)[C@H]1O)C(O)=O"  # (+)-quinic acid
print(is_quinic_acid(smiles_example))