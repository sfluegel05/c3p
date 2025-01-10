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
    
    # General long unsaturated carbon chain with some flexibility in double bond locations
    chain_pattern = Chem.MolFromSmarts("CCCCCCCCC(=C)C(=C)C(=C)CC")
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "Does not have recognizable long chain with multiple unsaturations"
    
    # Incorporate either an ethanolamine or glycerol-like structure
    ethanolamine_pattern = Chem.MolFromSmarts("NCCO")
    glycerol_pattern = Chem.MolFromSmarts("C(O)C(O)CO")
    if not (mol.HasSubstructMatch(ethanolamine_pattern) or mol.HasSubstructMatch(glycerol_pattern)):
        return False, "No characteristic endocannabinoid backbone found"
    
    # Look for either ester or amide linkage
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    if not (mol.HasSubstructMatch(ester_pattern) or mol.HasSubstructMatch(amide_pattern)):
        return False, "No ester or amide linkage found"
    
    # Check for additional features like hydroxyl or epoxy as optional
    hydroxyl_epoxy_pattern = Chem.MolFromSmarts("[OH] | [O]")
    if mol.HasSubstructMatch(hydroxyl_epoxy_pattern):
        return True, "Contains endocannabinoid-like features with hydroxyl/epoxy group"
    
    return True, "Contains endocannabinoid-like features"

# Test cases should include the SMILES strings provided above to verify the function works correctly.