"""
Classifies: CHEBI:67197 endocannabinoid
"""
from rdkit import Chem

def is_endocannabinoid(smiles: str):
    """
    Determines if a molecule is an endocannabinoid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an endocannabinoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Prepare more flexible patterns
    # Recognize long unsaturated carbon chains generally, allowing for flexibility
    chain_pattern = Chem.MolFromSmarts("C=C")
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "Does not have recognizable unsaturated carbon chain"
    
    # Recognize either ethanolamine or glycerol-like structure
    ethanolamine_pattern = Chem.MolFromSmarts("NCCO")
    glycerol_pattern = Chem.MolFromSmarts("C(O)C(O)C")
    if not (mol.HasSubstructMatch(ethanolamine_pattern) or mol.HasSubstructMatch(glycerol_pattern)):
        return False, "No characteristic endocannabinoid backbone found"

    # Look for either ester or amide linkage to ensure linkage to long chain
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    if not (mol.HasSubstructMatch(ester_pattern) or mol.HasSubstructMatch(amide_pattern)):
        return False, "No ester or amide linkage found"
    
    # Check for presence of hydroxyl or epoxy groups, optional but common in endocannabinoids
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    epoxy_pattern = Chem.MolFromSmarts("C1OC1")
    if mol.HasSubstructMatch(hydroxyl_pattern) or mol.HasSubstructMatch(epoxy_pattern):
        return True, "Contains endocannabinoid-like features with possible functional groups"

    return True, "Contains endocannabinoid-like features"

# The users can test the function with various endocannabinoid SMILES strings to verify accuracy.