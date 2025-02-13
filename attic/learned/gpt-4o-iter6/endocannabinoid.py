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

    # Check for presence of long carbon chain, with or without double bonds
    long_chain_pattern = Chem.MolFromSmarts("[C;!$(C(=O))]~[C;!$(C(=O))]~[C;!$(C(=O))]~[C;!$(C(=O))]")
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "No long carbon chain found typical of endocannabinoids."

    # Check for possible functional groups: ethanolamine, amide, or glycerol backbone
    ethanolamine_pattern = Chem.MolFromSmarts("NCCO")
    amide_pattern = Chem.MolFromSmarts("NCC(=O)")
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)CO")
    has_key_group = mol.HasSubstructMatch(ethanolamine_pattern) or mol.HasSubstructMatch(amide_pattern) or mol.HasSubstructMatch(glycerol_pattern)
    if not has_key_group:
        return False, "No ethanolamine, amide, or glycerol backbone found."

    # Check for ether, ester, or amide linkages
    ether_pattern = Chem.MolFromSmarts("OCC")
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    amide_linkage_pattern = Chem.MolFromSmarts("C(=O)N")
    if not (mol.HasSubstructMatch(ether_pattern) or mol.HasSubstructMatch(ester_pattern) or mol.HasSubstructMatch(amide_linkage_pattern)):
        return False, "No ether, ester, or amide linkages found."

    return True, "Matches characteristics of known endocannabinoids."