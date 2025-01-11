"""
Classifies: CHEBI:55493 1-O-acylglycerophosphoethanolamine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_1_O_acylglycerophosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-O-acylglycerophosphoethanolamine based on its SMILES string.

    A 1-O-acylglycerophosphoethanolamine is a glycerophosphoethanolamine having an O-acyl substituent at the 1-position of the glycerol fragment.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-O-acylglycerophosphoethanolamine, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for phosphoethanolamine group
    phosphoethanolamine_smarts = "O[P](O)(=O)OCCN"
    phosphoethanolamine_pattern = Chem.MolFromSmarts(phosphoethanolamine_smarts)
    if not mol.HasSubstructMatch(phosphoethanolamine_pattern):
        return False, "No phosphoethanolamine group found"

    # Look for ester groups (acyl chains attached via ester linkage)
    ester_pattern = Chem.MolFromSmarts("C(=O)O[C]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Expected one ester group, found {len(ester_matches)}"

    # Optionally, check for a free hydroxyl group (at position 2)
    hydroxyl_pattern = Chem.MolFromSmarts("[C;H1][OH]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 1:
        return False, "No free hydroxyl group at position 2 found"

    return True, "Contains phosphoethanolamine group and one acyl chain attached via ester linkage"

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': '1-O-acylglycerophosphoethanolamine',
        'definition': 'A glycerophosphoethanolamine having an unspecified O-acyl substituent at the 1-position of the glycerol fragment.',
        'parents': []
    }
}