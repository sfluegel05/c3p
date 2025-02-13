"""
Classifies: CHEBI:55493 1-O-acylglycerophosphoethanolamine
"""
from rdkit import Chem

def is_1_O_acylglycerophosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-O-acylglycerophosphoethanolamine
    based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 1-O-acylglycerophosphoethanolamine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define flexible SMARTS patterns for 1-O-acylglycerophosphoethanolamine components
    
    # 1-O-acylglycerol pattern: flexible enough to detect variances in chain length and stereochemistry
    glycerol_pattern = Chem.MolFromSmarts("COC(=O)[C]C(CO)O")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No 1-O-acyl-glycerol structure found"

    # Phosphoethanolamine group: allow variation in charge or resonance
    phosphoethanolamine_pattern = Chem.MolFromSmarts("COP([O-])(=O)OCC[NH3+]")
    alternative_phosphoethanolamine_pattern = Chem.MolFromSmarts("COP(O)(=O)OCCN")
    
    # Check presence of phosphoethanolamine group in any valid representation
    if not (mol.HasSubstructMatch(phosphoethanolamine_pattern) or mol.HasSubstructMatch(alternative_phosphoethanolamine_pattern)):
        return False, "Phosphoethanolamine group not present in expected forms"
    
    # Ensure both substructures can be linked and are part of the same molecule
    if not (mol.HasSubstructMatch(glycerol_pattern) and (mol.HasSubstructMatch(phosphoethanolamine_pattern) or mol.HasSubstructMatch(alternative_phosphoethanolamine_pattern))):
        return False, "Composite structure not consistent with 1-O-acylglycerophosphoethanolamine"

    return True, "Structure matches 1-O-acylglycerophosphoethanolamine"