"""
Classifies: CHEBI:17517 phosphatidylglycerol
"""
from rdkit import Chem

def is_phosphatidylglycerol(smiles: str):
    """
    Determines if a molecule is a phosphatidylglycerol based on its SMILES string.
    A phosphatidylglycerol has a glycerol backbone attached to a phosphate group and another glycerol moiety,
    often with two esterified fatty acid chains.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidylglycerol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for phosphate group linked to glycerol (P-O-C)
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)OCC(O)CO")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphatidylglycerol phosphate linkage found"
    
    # Check for esterified glycerol backbone (C(=O)OC)
    ester_pattern = Chem.MolFromSmarts("OC(=O)C")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester linkages, need at least 2"
    
    # Check for fatty acid chains (linear carbon chains)
    alkyl_chain_pattern = Chem.MolFromSmarts("C~C~C~C")
    alkyl_chain_matches = mol.GetSubstructMatches(alkyl_chain_pattern)
    if len(alkyl_chain_matches) < 2:
        return False, f"Insufficient alkyl chains for fatty acids, found {len(alkyl_chain_matches)}"

    return True, "Contains phosphatidylglycerol structural motifs"