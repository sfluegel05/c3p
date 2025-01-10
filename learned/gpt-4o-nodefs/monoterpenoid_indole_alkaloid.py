"""
Classifies: CHEBI:65323 monoterpenoid indole alkaloid
"""
from rdkit import Chem

def is_monoterpenoid_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid indole alkaloid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is a monoterpenoid indole alkaloid, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Enhanced SMARTS pattern for indole moiety including common variants.
    indole_patterns = [
        Chem.MolFromSmarts("c1c[nH]c2cccc3c2c1ccc3"),  # Standard indole
        Chem.MolFromSmarts("c1ccc2c(c1)[nH]c3c2cccc3"), # Substituted indole variants
    ]
    
    if not any(mol.HasSubstructMatch(pattern) for pattern in indole_patterns):
        return False, "No indole moiety found, or misidentified"

    # SMARTS pattern for a monoterpenoid moiety using isoprene-like units.
    isoprene_unit = Chem.MolFromSmarts("C(C)(C)C=C")

    # Check for presence of monoterpenoid features
    if not mol.HasSubstructMatch(isoprene_unit):
        return False, "No monoterpenoid features found"

    # Check for common monoterpenoid indole alkaloid functional groups like methoxy or esters
    methoxy_pattern = Chem.MolFromSmarts("[OX2]C")
    ester_pattern = Chem.MolFromSmarts("C(=O)OC")

    has_methoxy = mol.HasSubstructMatch(methoxy_pattern)
    has_ester = mol.HasSubstructMatch(ester_pattern)

    if not (has_methoxy or has_ester):
        return False, "Lacks common functional groups of monoterpenoid indole alkaloids"

    return True, "Contains indole moiety and monoterpenoid features with relevant functional groups"