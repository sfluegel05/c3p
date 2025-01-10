"""
Classifies: CHEBI:26004 phenylpropanoid
"""
from rdkit import Chem

def is_phenylpropanoid(smiles: str):
    """
    Determines if a molecule is a phenylpropanoid based on its SMILES string.
    Phenylpropanoids are organic aromatic compounds with a phenylpropane skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phenylpropanoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define phenylpropanoid pattern (aromatic phenyl ring with three-carbon propenyl group -C3H5)
    phenylpropanoid_pattern = Chem.MolFromSmarts("c1ccccc1CCC")
    if mol.HasSubstructMatch(phenylpropanoid_pattern):
        return True, "Contains phenylpropane skeleton"
    
    # Additional checks for known phenylpropanoid subclasses could be added here
    # such as checking for flavonoids, coumarins, etc.

    return False, "Does not match phenylpropanoid structure"