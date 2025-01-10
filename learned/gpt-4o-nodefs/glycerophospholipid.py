"""
Classifies: CHEBI:37739 glycerophospholipid
"""
from rdkit import Chem

def is_glycerophospholipid(smiles: str):
    """
    Determines if a molecule is a glycerophospholipid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycerophospholipid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone pattern (C-C-C with 2 or more oxygens attached)
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Two long acyl chains (usually represented by -C(=O)CCC-), so a general carbon chain
    acyl_chain_pattern = Chem.MolFromSmarts("C(=O)CC")
    acyl_matches = mol.GetSubstructMatches(acyl_chain_pattern)
    if len(acyl_matches) < 2:
        return False, "Less than two acyl chains found"

    # Phosphate group pattern (P(=O)(O)O[R])
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)O")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"

    return True, "Contains glycerol backbone with two acyl chains and a phosphate group"

# Additional logic and checking could be added here to further confirm the identity based on subclass