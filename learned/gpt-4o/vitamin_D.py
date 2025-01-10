"""
Classifies: CHEBI:27300 vitamin D
"""
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcNumRotatableBonds
from rdkit.Chem.rdmolops import SanitizeFlags

def is_vitamin_D(smiles: str):
    """
    Determines if a molecule is a vitamin D compound based on its SMILES string.
    Includes tests for open B-ring and typical secosteroid configurations.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a vitamin D compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    try:
        Chem.SanitizeMol(mol, sanitizeOps=SanitizeFlags.SANITIZE_ALL)

        # Adjusted SMARTS pattern for open B-ring with variations
        open_b_ring_pattern = Chem.MolFromSmarts("[C@H]1(C=C/C2CCC[C@]3(C=C/C=C/4[C@]5([C@@]2(CCC3)C)CC(CCC5)CC4=C)[H])[H]")
        if not mol.HasSubstructMatch(open_b_ring_pattern):
            return False, "No or atypical secosteroid-like open B-ring structure found"
        
        # Identify conjugated triene that permits some deviations
        triene_pattern = Chem.MolFromSmarts("C=C/C=C/C=C")
        if not mol.HasSubstructMatch(triene_pattern):
            return False, "No conjugated triene system detected"
        
        # Detect at least two or more hydroxyl groups, typical yet variable in vitamin D
        hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
        hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
        if len(hydroxyl_matches) < 2:
            return False, "Less than two hydroxyl groups detected"

        # Additional checks for secosteroid features might include flexibility in structure
        # Based on rotatable bonds (flexible feature check) assumed in secosteroid structures
        if CalcNumRotatableBonds(mol) < 3:
            return False, "Insufficient structural flexibility for vitamin D"

    except Exception as e:
        return False, f"SMILES parsing or feature detection failed: {str(e)}"

    # If all checks are passed, classify as vitamin D
    return True, "Matches typical vitamin D structural features"