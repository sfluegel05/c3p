"""
Classifies: CHEBI:26493 quinic acid
"""
from rdkit import Chem

def is_quinic_acid(smiles: str):
    """
    Determines if a molecule is a quinic acid or its derivative based on its SMILES string.
    Quinic acid is a cyclitol carboxylic acid, often found in its esterified form.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a quinic acid or relevant derivative, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for a cyclohexane ring
    ring_info = mol.GetRingInfo()
    if not any(len(ring) == 6 for ring in ring_info.AtomRings()):
        return False, "No cyclohexane ring found"

    # Check hydroxy groups, focusing on stereochemistry specific to quinic acid:
    hydroxy_pattern = Chem.MolFromSmarts("[C@@H](O)")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    if len(hydroxy_matches) < 3:
        return False, f"Not enough hydroxyl groups with correct stereochemistry found, got {len(hydroxy_matches)}"

    # Check for carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Check for ester modifications, common in quinic derivatives:
    ester_pattern = Chem.MolFromSmarts("OC=O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) > 0:
        return True, "Contains features consistent with esterified quinic acid"

    return True, "Contains features consistent with quinic acid"

# Example use case
smiles_example = "O[C@H]1C[C@@](O)(C[C@H](O)[C@H]1O)C(O)=O"  # (+)-quinic acid
print(is_quinic_acid(smiles_example))