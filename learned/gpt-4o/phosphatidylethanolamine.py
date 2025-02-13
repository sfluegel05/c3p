"""
Classifies: CHEBI:16038 phosphatidylethanolamine
"""
from rdkit import Chem

def is_phosphatidylethanolamine(smiles: str):
    """
    Determines if a molecule is a phosphatidylethanolamine based on its SMILES string.
    A phosphatidylethanolamine is a glycerophospholipid where two fatty acids are esterified to a glycerol backbone,
    and a phosphatidyl group is esterified to ethanolamine.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidylethanolamine, False otherwise.
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone
    glycerol_pattern = Chem.MolFromSmarts("[C@H](COP(O)(=O)OCCN)O")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol-like backbone with phosphoethanolamine group attached found"
        
    # More flexible ester linkage patterns (to include ester diversity)
    ester_pattern = Chem.MolFromSmarts("COC(=O)")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester linkages, need at least 2"

    # Check for phosphoethanolamine group, allowing for N-methyl groups
    phosphoethanolamine_pattern = Chem.MolFromSmarts("P(OCC[NH2,NH3+,NH(CH3),NH(CH3)2])")
    if not mol.HasSubstructMatch(phosphoethanolamine_pattern):
        return False, "No phosphoethanolamine or derived group found"
    
    # Successfully matched all criteria
    return True, "Contains glycerol backbone with 2 esterified fatty acids and a phosphoethanolamine group"