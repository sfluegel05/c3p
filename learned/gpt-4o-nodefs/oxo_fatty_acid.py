"""
Classifies: CHEBI:59644 oxo fatty acid
"""
from rdkit import Chem

def is_oxo_fatty_acid(smiles: str):
    """
    Determines if a molecule is an oxo fatty acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an oxo fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxylic acid group (-C(=O)O)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Look for additional ketone group (not part of carboxylic acid)
    ketone_pattern = Chem.MolFromSmarts("[CX3](=O)[C;!$(C(=O)O)]")
    if not mol.HasSubstructMatch(ketone_pattern):
        return False, "No additional ketone group found"

    # Check for a bigger carbon skeleton; allow for more carbons, rings, and unsaturations
    # This pattern now allows for any arrangement of carbons while ensuring a long chain backbone
    carbon_chain_pattern = Chem.MolFromSmarts("C-C-C-C-C-C-C")  # flexible pattern allowing more configurations
    if not mol.HasSubstructMatch(carbon_chain_pattern):
        return False, "Insufficient carbon backbone"

    # Optionally verify unsaturations or cyclic components in more complex molecules
    if Chem.MolFromSmarts("C=C") and mol.HasSubstructMatch(carboxylic_acid_pattern):
        unsaturation_pattern = Chem.MolFromSmarts("C=C")
        if not mol.HasSubstructMatch(unsaturation_pattern):
            return False, "Lacks unsaturation often present in oxo fatty acids"
     
    return True, "Contains carboxylic acid and additional ketone group with sufficient carbon backbone"