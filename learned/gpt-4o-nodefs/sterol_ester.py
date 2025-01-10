"""
Classifies: CHEBI:35915 sterol ester
"""
from rdkit import Chem

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

    # Sterol core identification: typical core structure involving four fused rings
    # Base steroid nucleus, allowing for substitution and stereochemistry variability
    sterol_smarts = "[#6]12[#6][#6]3[#6]([#6]1)[#6][#6]4[C@]2([#6])CC[C@]45C[C@H]3CC[C@@H]5C"
    sterol_core = Chem.MolFromSmarts(sterol_smarts)

    if not mol.HasSubstructMatch(sterol_core):
        return False, "No sterol core structure found"
        
    # Ester group pattern (-C(=O)O-)
    ester_pattern_smarts = "[#6](=O)O"  # Valid for ester linkage
    ester_pattern = Chem.MolFromSmarts(ester_pattern_smarts)
    
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester linkage found"

    # Typically, sterol esters have long chains; count carbons in chains attached to ester
    aliphatic_chain_length = 10  # Arbitrary length suggesting a fatty acid
    chain_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if len(chain_atoms) < aliphatic_chain_length:
        return False, "Insufficient aliphatic chain found"

    return True, "Contains sterol core structure esterified with a sufficient long aliphatic chain"

# Example metadata
__metadata__ = {   'chemical_class': {   'id': 'CHEBI:None',
                          'name': 'sterol ester',
                          'definition': 'Chemical entities consisting of a sterol esterified with a fatty acid'},
    'config': {   'llm_model_name': 'rdkit_based_model',
                  'f1_threshold': 0.8,
                  'max_attempts': 5},
    'message': None,
    'attempt': 2,
    'success': True}