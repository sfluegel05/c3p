"""
Classifies: CHEBI:35346 11beta-hydroxy steroid
"""
from rdkit import Chem

def is_11beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is an 11beta-hydroxy steroid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an 11beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic steroid backbone pattern: 17 carbon atoms arranged in four fused rings
    # This example uses a general steroid core scaffold (cyclopenta[a]phenanthrene) as a simple pattern
    steroid_pattern = Chem.MolFromSmarts('C1CCC2C(C1)CCC3C2CCC4C3CCC4') 
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "Basic steroid backbone not found"

    # Look for the 11beta-hydroxy group: CCOH at position 11 with beta configuration
    # Simplified SMARTS assuming the general placement
    # Beta configuration commonly signifies a distinct stereochemistry
    hydroxy_11beta_pattern = Chem.MolFromSmarts('[C@@H](O)[C@]1CCCCC1')
    if not mol.HasSubstructMatch(hydroxy_11beta_pattern):
        return False, "No 11beta-hydroxy group found with correct configuration"

    return True, "Contains basic steroid backbone with an 11beta-hydroxy group of the correct configuration"