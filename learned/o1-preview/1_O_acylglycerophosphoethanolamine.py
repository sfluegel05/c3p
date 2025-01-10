"""
Classifies: CHEBI:55493 1-O-acylglycerophosphoethanolamine
"""
from rdkit import Chem

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

    # Define SMARTS for phosphoethanolamine group (allowing for different protonation states)
    phosphoethanolamine_smarts = 'O[P](=O)(O)-O-CC[N]'
    phosphoethanolamine = Chem.MolFromSmarts(phosphoethanolamine_smarts)
    if not mol.HasSubstructMatch(phosphoethanolamine):
        return False, "Phosphoethanolamine group not found"

    # Define SMARTS for glycerol backbone with free hydroxyl at position 2
    glycerol_smarts = 'OCC(O)CO'
    glycerol = Chem.MolFromSmarts(glycerol_smarts)
    if not mol.HasSubstructMatch(glycerol):
        return False, "Glycerol backbone with free hydroxyl at position 2 not found"

    # Define SMARTS for ester linkage at position 1 (acyl group)
    ester_smarts = 'OC(=O)[C]'
    ester = Chem.MolFromSmarts(ester_smarts)
    if not mol.HasSubstructMatch(ester):
        return False, "Ester linkage at position 1 not found"

    # Ensure correct connectivity by defining a full pattern
    full_pattern_smarts = 'OC(=O)C(O)COP(O)(=O)OCC[N]'
    full_pattern = Chem.MolFromSmarts(full_pattern_smarts)
    if not mol.HasSubstructMatch(full_pattern):
        return False, "Molecule does not match the complete 1-O-acylglycerophosphoethanolamine structure"

    return True, "Matches 1-O-acylglycerophosphoethanolamine structure"