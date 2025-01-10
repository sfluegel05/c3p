"""
Classifies: CHEBI:17962 CDP-diacylglycerol
"""
from rdkit import Chem

def is_CDP_diacylglycerol(smiles: str):
    """
    Determines if a molecule is a CDP-diacylglycerol based on its SMILES string.
    A CDP-diacylglycerol is characterized by a glycerol backbone with two acyl chains 
    at positions 1 and 2, and a cytidine diphosphate group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a CDP-diacylglycerol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generalized pattern for the CDP group
    # This pattern searches for cytidine with two phosphates, allowing variations
    # using flexible phospho linkages and positional adaptations
    cdp_smarts = "n1c(COC(OP(O)=O)P(O)=O)[C@H]([C@@H]1O)COP(O)(O)=O"
    cdp_pattern = Chem.MolFromSmarts(cdp_smarts)
    if not mol.HasSubstructMatch(cdp_pattern):
        return False, "No CDP (Cytidine diphosphate) group found"

    # Look for glycerol backbone allowing for some stereochemical variability
    glycerol_pattern = Chem.MolFromSmarts("[C@H](O)C(CO)O")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Look for two ester groups (-C(=O)O-), specifically at glycerol 1,2 positions.
    ester_pattern = Chem.MolFromSmarts("C(=O)O(C(O)CO)")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} acyl chains, need at least 2"

    # Verify that ester patterns are connected to glycerol correct positions
    attachment_positions = {match[1] for match in ester_matches}
    if len(attachment_positions) < 2:
        return False, "Acyl chains not correctly placed at 1,2-glycerol positions"

    # Verify the length criteria of acyl chains to ensure they are sufficiently long
    acyl_chain_count = sum(
        1 for atom in mol.GetAtoms()
        if atom.GetAtomicNum() == 6 and len(atom.GetBonds()) in [2, 3]
    )
    
    if acyl_chain_count < 20:
        return False, "Acyl chains too short to classify"

    return True, "Contains CDP group, appropriate glycerol backbone, and properly attached acyl chains"