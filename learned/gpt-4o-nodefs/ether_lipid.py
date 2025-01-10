"""
Classifies: CHEBI:64611 ether lipid
"""
from rdkit import Chem

def is_ether_lipid(smiles: str):
    """
    Determines if a molecule is an ether lipid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an ether lipid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Adjusted pattern for ether linkage
    ether_pattern = Chem.MolFromSmarts("COC")
    if not mol.HasSubstructMatch(ether_pattern):
        return False, "No ether linkage found"

    # More flexible glycerol backbone patterns
    glycerol_patterns = [
        Chem.MolFromSmarts("C(CO)CO"),  # Traditional glycerol pattern
        Chem.MolFromSmarts("C(C)CO"),   # Slightly modified for ether lipids
        Chem.MolFromSmarts("OCC(O)C")   # Allows for different orientations
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in glycerol_patterns):
        return False, "No suitable glycerol backbone found"
    
    # Recognize polar head groups (consider common ones)
    head_group_patterns = [
        Chem.MolFromSmarts("OP(O)(=O)O"),  # Phosphate group
        Chem.MolFromSmarts("N(C)(C)C")    # Choline headgroup (common in phosphocholines)
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in head_group_patterns):
        return False, "No recognized polar head group found (e.g., phosphate, choline)"

    # Verify presence of long hydrocarbon chains
    hydrocarbon_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]") 
    if not mol.HasSubstructMatch(hydrocarbon_chain_pattern):
        return False, "No adequate hydrocarbon chain found"

    return True, "Identified ether linkage, glycerol backbone, polar head group, and hydrocarbon chains typical of ether lipids"