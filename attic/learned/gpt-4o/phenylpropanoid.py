"""
Classifies: CHEBI:26004 phenylpropanoid
"""
from rdkit import Chem

def is_phenylpropanoid(smiles: str):
    """
    Determines if a molecule is a phenylpropanoid based on its SMILES string.
    Phenylpropanoids are organic aromatic compounds with a phenylpropane-derived skeleton.

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

    # Create SMARTS patterns that represent common features in phenylpropanoids
    phenylpropanoid_pattern1 = Chem.MolFromSmarts("c1ccccc1CC")  # Example: simple phenylpropane
    flavonoid_pattern = Chem.MolFromSmarts("c1c(=O)oc2ccccc2c1") # Example: flavone
    coumarin_pattern = Chem.MolFromSmarts("c1ccc2oc1ccc2=O")  # Example: coumarin

    # Check for phenylpropanoid-like patterns
    if mol.HasSubstructMatch(phenylpropanoid_pattern1):
        return True, "Contains phenylpropane skeleton"
    if mol.HasSubstructMatch(flavonoid_pattern):
        return True, "Contains flavone-like structure"
    if mol.HasSubstructMatch(coumarin_pattern):
        return True, "Contains coumarin-like structure"

    # Additional logic for derivatives or structurally similar compounds,
    # which can improve hit rate over simply phenylpropane and can include
    # naturally occurring substituents, fused aromatic systems, etc.

    return False, "Does not match phenylpropanoid structure"