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

    # Check for ether linkage pattern
    ether_pattern = Chem.MolFromSmarts("COC")
    if not mol.HasSubstructMatch(ether_pattern):
        return False, "No ether linkage found"

    # Check for glycerol backbone
    glycerol_pattern = Chem.MolFromSmarts("C(C)(CO)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Identify phosphate groups
    phosphate_pattern = Chem.MolFromSmarts("OP(O)(=O)O")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"

    # Check for long hydrocarbon chains (at least one)
    # Here we match chains with at least 8 carbons as a minimal length threshold.
    hydrocarbon_chain_pattern = Chem.MolFromSmarts("[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]") 
    if not mol.HasSubstructMatch(hydrocarbon_chain_pattern):
        return False, "No long hydrocarbon chain found"
    
    return True, "Contains ether linkage with glycerol backbone, phosphate group, and long hydrocarbon chains"

__metadata__ = {   
    'chemical_class': {   
        'name': 'ether lipid',
        'definition': 'A lipid containing an ether bond, a glycerol backbone, '
                      'and typically a phosphate group with long hydrocarbon chains.'
    }
}