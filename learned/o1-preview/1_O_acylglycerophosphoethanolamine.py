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

    # Define SMARTS patterns

    # Phosphoethanolamine group at position 3
    phosphoethanolamine_smarts = "O[P](=O)(O)OCCN"
    phosphoethanolamine = Chem.MolFromSmarts(phosphoethanolamine_smarts)
    if phosphoethanolamine is None:
        return False, "Invalid phosphoethanolamine SMARTS pattern"

    # Check for phosphoethanolamine group
    if not mol.HasSubstructMatch(phosphoethanolamine):
        return False, "Phosphoethanolamine group not found at position 3"

    # Ester linkage at position 1 of glycerol backbone
    ester_smarts = "OC(=O)[C]"
    ester = Chem.MolFromSmarts(ester_smarts)
    if ester is None:
        return False, "Invalid ester SMARTS pattern"

    # Check for ester linkage
    if not mol.HasSubstructMatch(ester):
        return False, "No ester linkage at position 1 found"

    # Glycerol backbone with free hydroxyl at position 2
    glycerol_smarts = "OCC(O)CO"
    glycerol = Chem.MolFromSmarts(glycerol_smarts)
    if glycerol is None:
        return False, "Invalid glycerol SMARTS pattern"

    # Check for glycerol backbone
    if not mol.HasSubstructMatch(glycerol):
        return False, "Glycerol backbone with free hydroxyl at position 2 not found"

    # Ensure that the ester linkage, glycerol backbone, and phosphoethanolamine group are connected correctly
    # Combine patterns to check proper connectivity
    full_pattern_smarts = "OC(=O)[C@@H]CO[P](=O)(O)OCCN"
    full_pattern = Chem.MolFromSmarts(full_pattern_smarts)
    if full_pattern is None:
        return False, "Invalid full SMARTS pattern"

    if not mol.HasSubstructMatch(full_pattern):
        # Try without stereochemistry
        full_pattern_smarts = "OC(=O)C(O)COP(=O)(O)OCCN"
        full_pattern = Chem.MolFromSmarts(full_pattern_smarts)
        if full_pattern is None:
            return False, "Invalid full SMARTS pattern without stereochemistry"
        if not mol.HasSubstructMatch(full_pattern):
            return False, "Molecule does not match the complete 1-O-acylglycerophosphoethanolamine structure"

    return True, "Matches 1-O-acylglycerophosphoethanolamine structure"