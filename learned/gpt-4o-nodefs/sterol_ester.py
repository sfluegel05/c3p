"""
Classifies: CHEBI:35915 sterol ester
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sterol_ester(smiles: str):
    """
    Determines if a molecule is a sterol ester based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sterol ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Refined sterol core pattern (steroid nucleus with characteristic rings and stereochemistry)
    sterol_smarts = "[#6]1[#6][#6]2[#6]3[#6]([#6]1)[#6][#6][#6]4[#6]([#6]3)[#6][#6][#6]5[#6](C2)[#6](C([#6][#6]5)[C@@H]4)O"
    sterol_core = Chem.MolFromSmarts(sterol_smarts)
    
    if not mol.HasSubstructMatch(sterol_core):
        return False, "No sterol core structure found"
        
    # Ester group pattern (-C(=O)O-)
    ester_pattern_smarts = "C(=O)O"
    ester_pattern = Chem.MolFromSmarts(ester_pattern_smarts)
    
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester linkage found"

    # Check for long aliphatic chain (using flexible chain detection)
    carbon_chain_smarts = "[#6]~[#6]~[#6]~[#6]~[#6]~[#6]" # Arbitrary extension for long chain
    carbon_chain = Chem.MolFromSmarts(carbon_chain_smarts)
    
    if not mol.HasSubstructMatch(carbon_chain):
        return False, "No sufficient aliphatic chain found"

    return True, "Contains sterol core structure esterified with a long aliphatic chain"

# Example metadata
__metadata__ = {   'chemical_class': {   'id': 'CHEBI:None',
                          'name': 'sterol ester',
                          'definition': 'Chemical entities consisting of a sterol esterified with a fatty acid'},
    'config': {   'llm_model_name': 'rdkit_based_model',
                  'f1_threshold': 0.8,
                  'max_attempts': 5},
    'message': None,
    'attempt': 1,
    'success': True}