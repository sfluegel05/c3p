"""
Classifies: CHEBI:35343 17beta-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_17beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 17beta-hydroxy steroid based on its SMILES string.
    A 17beta-hydroxy steroid is defined by having a steroid backbone with a beta-configured hydroxy group at position 17.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 17beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # General steroid backbone pattern - four rings A/B/C/D with flexibility in stereo and connectivity
    # These SMARTS patterns aim to capture the core steroid backbone, allowing for stereochemical variations
    # The specific SMARTS patterns should be developed in accordance with the structure of known steroids
    steroid_patterns = [
        Chem.MolFromSmarts("C1C2CC3CC4=C(C=CC=C4)CC3C2C1"),  # Simplified steroid backbone
        Chem.MolFromSmarts("C1C2CC3CC4=CC=CC=C4CC3C2C1"),   # A more generic pattern capturing aromatic C-ring possibility
    ]
    
    if not any(mol.HasSubstructMatch(pattern) for pattern in steroid_patterns):
        return False, "No suitable steroid backbone found"
    
    # Look for the 17beta-hydroxy group configuration specifically
    # Assuming that the 17th carbon is part of the general steroid alkane ring system,
    # Adjust the pattern to accurately represent the possible locations and include stereochemistry.
    # This SMARTS assumes we know the general 17-OH placement chemistry typical in steroids.
    hydroxy_pattern = Chem.MolFromSmarts("[C@@H](O)")  # An stereospecific pattern assuming beta configuration
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "No 17beta-hydroxy group detected"
    
    return True, "Detected a 17beta-hydroxy steroid configuration"