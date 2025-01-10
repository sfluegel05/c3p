"""
Classifies: CHEBI:27300 vitamin D
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_vitamin_D(smiles: str):
    """
    Determines if a molecule is a vitamin D compound based on its SMILES string.
    Vitamin D compounds are defined as fat-soluble hydroxy seco-steroids.

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
    
    # Define an improved seco-steroid pattern (as a vitamin D specific structure)
    # The core seco-steroid vitamin D structure with a broken ring
    # We will simplify the matching of a typical vitamin D seco structure
    seco_steroid_pattern = Chem.MolFromSmarts("C1CCC2C=C(C)CCC2C1")  # Simpler pattern matching B-ring

    if seco_steroid_pattern is None:
        return (None, None)  # If the pattern itself is not correctly constructed, prevent failure
    
    if not mol.HasSubstructMatch(seco_steroid_pattern):
        return False, "No seco-steroid backbone typical of vitamin D identified"
    
    # Check for presence of hydroxyl groups (-OH)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")  # OH group
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 2:
        return False, f"Only {len(hydroxyl_matches)} hydroxyl groups found, typical vitamin D compounds need more"

    # Check molecular properties indicative of lipophilicity
    mol_logP = rdMolDescriptors.CalcCrippenDescriptors(mol)[0]
    if mol_logP < 3:  # Tune threshold if needed for specific definition
        return False, "Molecule not hydrophobic enough for vitamin D classification"
   
    return True, "Matches typical vitamin D seco-steroid structure with hydroxyl groups"