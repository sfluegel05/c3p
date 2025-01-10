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

    # Long unsaturated carbon chain pattern
    # This pattern looks for long carbon chains with multiple unsaturations
    long_chain_pattern = Chem.MolFromSmarts("C=C-@C=C-@C=C(@C=C-@C=C)CCCC")
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "Does not have long unsaturated carbon chains"

    # Ethanoloamine group pattern -NCCO
    ethanolamine_pattern = Chem.MolFromSmarts("NCCO")
    if not mol.HasSubstructMatch(ethanolamine_pattern):
        return False, "No ethanolamine group found"
    
    # Glycerol backbone pattern
    glycerol_pattern = Chem.MolFromSmarts("C(O)C(O)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol or similar backbone found"

    # Ester or amide linkage pattern
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    if not (mol.HasSubstructMatch(ester_pattern) or mol.HasSubstructMatch(amide_pattern)):
        return False, "No ester or amide linkage found"
    
    # Optional: check for hydroxyl or epoxy groups (not strictly necessary for all endocannabinoids)
    hydroxyl_epoxy_pattern = Chem.MolFromSmarts("[OH,O]")
    if mol.HasSubstructMatch(hydroxyl_epoxy_pattern):
        return True, "Contains endocannabinoid-like features and hydroxyl/epoxy group"

    return True, "Contains endocannabinoid-like features"