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
        # Sanitize molecule (previously applied but omitted; RDKit usually autoinfers sanitization need)
        Chem.SanitizeMol(mol)
        
        # Corrected SMARTS pattern for typical secosteroid-like structure
        open_b_ring_pattern = Chem.MolFromSmarts("C1C=C2CCCCC2=C1")
        if not mol.HasSubstructMatch(open_b_ring_pattern):
            return False, "No typical secosteroid-like open B-ring structure found"
        
        # Identify conjugated triene system with flexibility
        triene_pattern = Chem.MolFromSmarts("C=C-C=C-C=C")
        if not mol.HasSubstructMatch(triene_pattern):
            return False, "No conjugated triene system detected"
        
        # Detect two or more hydroxyl groups
        hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
        hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
        if len(hydroxyl_matches) < 2:
            return False, "Less than two hydroxyl groups detected"

        # Structural flexibility in secosteroid structure characterized by sufficient rotatable bonds
        if CalcNumRotatableBonds(mol) < 3:
            return False, "Insufficient structural flexibility for vitamin D"

    except Exception as e:
        return False, f"SMILES parsing or feature detection failed: {str(e)}"

    # If all checks are passed, classify as vitamin D
    return True, "Matches typical vitamin D structural features"