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
    
    # Potential Substructure Matches for Endocannabinoids

    # Longer polyunsaturated carbon chain
    extended_chain_pattern = Chem.MolFromSmarts("C=C.C=C.C=C")
    if not mol.HasSubstructMatch(extended_chain_pattern):
        return False, "Does not have recognizable polyunsaturated carbon chain"
    
    # Backbone patterns: ethanolamine or glycerol-like linkage
    ethanolamine_pattern = Chem.MolFromSmarts("NCCO")
    glycerol_pattern = Chem.MolFromSmarts("C(O)C(O)C")
    if not (mol.HasSubstructMatch(ethanolamine_pattern) or mol.HasSubstructMatch(glycerol_pattern)):
        return False, "No characteristic endocannabinoid backbone found"
    
    # Attachment via ester or amide linkage
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    if not (mol.HasSubstructMatch(ester_pattern) or mol.HasSubstructMatch(amide_pattern)):
        return False, "No ester or amide linkage found"
    
    # Optional but common functional groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    epoxy_pattern = Chem.MolFromSmarts("C1OC1")
    if mol.HasSubstructMatch(hydroxyl_pattern) or mol.HasSubstructMatch(epoxy_pattern):
        return True, "Contains endocannabinoid-like features with possible functional groups"

    return True, "Contains endocannabinoid-like features"

# The users can test the function with various endocannabinoid SMILES strings to verify accuracy.