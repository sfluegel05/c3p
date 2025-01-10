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

    # Accurate pattern for the CDP (cytidine diphosphate) group
    cdp_smarts = "n1c(COP(O)(=O)OP(O)(=O)O)cnc1C2OC(CO)C(O)C2O"
    cdp_pattern = Chem.MolFromSmarts(cdp_smarts)
    if not mol.HasSubstructMatch(cdp_pattern):
        return False, "No CDP (Cytidine diphosphate) group found"

    # Flexible glycerol backbone pattern
    glycerol_pattern = Chem.MolFromSmarts("OCC(CO)(O)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Look for two ester groups (-C(=O)O-)
    ester_pattern = Chem.MolFromSmarts("C(=O)OC[C@H](CO)O")  # Match positions in glycerol
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 2"

    # Check carbon chain length in acyl groups
    carbon_chain_length = sum(
        1 for atom in mol.GetAtoms()
        if atom.GetAtomicNum() == 6 and len(atom.GetNeighbors()) > 1
    )
    if carbon_chain_length < 20:
        return False, f"Acyl chains too short, carbon count {carbon_chain_length}"

    # Verify unique attachment positions in glycerol
    attachment_positions = {match[1] for match in ester_matches}
    if len(attachment_positions) < 2:
        return False, "Acyl chains not correctly attached at glycerol 1,2 positions"

    return True, "Contains CDP structure, glycerol backbone, and suitable acyl chains"