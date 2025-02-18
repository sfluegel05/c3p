"""
Classifies: CHEBI:50523 butenolide
"""
"""
Classifies: CHEBI:35489 butenolide
A gamma-lactone that consists of a 2-furanone skeleton and its substituted derivatives.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_butenolide(smiles: str):
    """
    Determines if a molecule is a butenolide based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a butenolide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for 2-furanone skeleton
    furanone_pattern = Chem.MolFromSmarts("O=C1OC=CC1")
    if not mol.HasSubstructMatch(furanone_pattern):
        return False, "No 2-furanone skeleton found"
    
    # Check for substitutions on the furanone ring
    substituted_furanone_pattern = Chem.MolFromSmarts("O=C1OC=C[C@@]1[!#1]")
    if not mol.HasSubstructMatch(substituted_furanone_pattern):
        return False, "2-furanone ring is not substituted"
    
    # Check for specific substitution patterns
    alkyl_chain_pattern = Chem.MolFromSmarts("[CX4][CX4][CX4][CX4]")
    alkoxy_pattern = Chem.MolFromSmarts("[OX2][CX4]")
    ester_pattern = Chem.MolFromSmarts("[OX2]C(=O)[CX4]")
    
    alkyl_chain_matches = mol.GetSubstructMatches(alkyl_chain_pattern)
    alkoxy_matches = mol.GetSubstructMatches(alkoxy_pattern)
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    substitution_patterns = []
    if alkyl_chain_matches:
        substitution_patterns.append("long alkyl chain")
    if alkoxy_matches:
        substitution_patterns.append("alkoxy group")
    if ester_matches:
        substitution_patterns.append("ester group")
    
    if not substitution_patterns:
        return False, "No characteristic butenolide substitution patterns found"
    
    reason = "Contains 2-furanone skeleton with the following substitutions: " + ", ".join(substitution_patterns)
    return True, reason