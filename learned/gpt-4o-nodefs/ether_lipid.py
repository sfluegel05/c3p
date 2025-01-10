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

    # Detect ether linkages - flexible atom types but with carbon chains
    ether_linkage_pattern = Chem.MolFromSmarts("[C]O[C]")  # Alkyl-O-Alkyl pattern
    if not mol.HasSubstructMatch(ether_linkage_pattern):
        return False, "No ether linkage pattern found"
    
    # Detect glycerol backbone with option for ether or ester linkages
    glycerol_patterns = [
        Chem.MolFromSmarts("[C@@H](CO[C])O[C]"),  # sn-glycero-3-phosphate form
        Chem.MolFromSmarts("[C@H](CO[C])O[C]"),   # alternative stereo
        Chem.MolFromSmarts("CO[C@H](CO)CO"),      # alternative etherified glycerol
        Chem.MolFromSmarts("[C@H](O[C])CO[C]")    # generalized form
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in glycerol_patterns):
        return False, "No identifiable glycerol backbone found"
    
    # Broad recognizing patterns for polar head groups
    head_group_patterns = [
        Chem.MolFromSmarts("OP(=O)(O)O"),  # Phosphate group
        Chem.MolFromSmarts("N(C)(C)C"),    # Choline
        Chem.MolFromSmarts("OC(=O)O"),     # Carboxyl group
        Chem.MolFromSmarts("P(=O)(O)OC"),  # Phosphonate group
        Chem.MolFromSmarts("NP(=O)(O)OC")  # Aminophosphate pattern
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in head_group_patterns):
        return False, "No recognizable polar head group found (e.g., phosphate, choline, carboxyl)"
    
    # Ensure long hydrocarbon chains, potentially unsaturated
    long_chain_pattern = Chem.MolFromSmarts("C(C)C(C)C")  # General repeated -CH2- units
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "No adequate long hydrocarbon chains found typical of ether lipids"
    
    return True, "Identified ether linkage, glycerol backbone, polar head group, and hydrocarbon chains typical of ether lipids"