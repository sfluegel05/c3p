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

    # Specific ether linkage pattern
    ether_linkage_pattern = Chem.MolFromSmarts("[!#1]-O-[!#1]")  # Non-hydrogen atoms on both sides of O
    if not mol.HasSubstructMatch(ether_linkage_pattern):
        return False, "No ether linkage pattern found"
    
    # Broad range of possible glycerol backbone patterns
    glycerol_patterns = [
        Chem.MolFromSmarts("C(CO)CO"),  # Standard glycerol
        Chem.MolFromSmarts("[C@@H](CO)CO"),  # Stereospecific
        Chem.MolFromSmarts("[C@H](CO)CO")   # Alternative stereo
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in glycerol_patterns):
        return False, "No suitable glycerol backbone found"
    
    # Expanded recognition of polar head groups
    head_group_patterns = [
        Chem.MolFromSmarts("OP(=O)(O)O"),  # Phosphate group
        Chem.MolFromSmarts("N(C)(C)C"),    # Choline
        Chem.MolFromSmarts("OC(=O)O"),     # carboxyl group
        Chem.MolFromSmarts("NP([O-])(=O)O") # phosphate - ammonia pattern
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in head_group_patterns):
        return False, "No recognized polar head group found (e.g., phosphate, choline, carboxyl)"
    
    # Recognition of long alkyl chains, possibly involving unsaturation
    long_chain_pattern = Chem.MolFromSmarts("[C](~[C,!H0])~[!#1]")  # General long carbon chain with possible unsaturation
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "No adequate long hydrocarbon chains found typical of ether lipids"
    
    return True, "Identified ether linkage, glycerol backbone, polar head group, and hydrocarbon chains typical of ether lipids"