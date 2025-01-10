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
    # This pattern searches for cytidine with two phosphates
    # The phosphate groups must be connected to cytidine in a typical way
    cdp_smarts = "n1c(COP(=O)(O)OP(=O)(O)O)[C@H]([C@@H]1O)CO"  # Adjusted general recognition pattern
    cdp_pattern = Chem.MolFromSmarts(cdp_smarts)
    if not mol.HasSubstructMatch(cdp_pattern):
        return False, "No CDP (Cytidine diphosphate) group found"

    # Look for glycerol backbone allowing for some stereochemical variability
    glycerol_pattern = Chem.MolFromSmarts("C(CO)(O)CO")  # Simplified for backbone detection
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Look for two acyl chain attachments (-C(=O)O), should be specifically at glycerol 1,2 positions.
    ester_pattern = Chem.MolFromSmarts("C(=O)OCC(CO)O")  # Adjusted to capture participation in ester linkages
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} acyl chains, need at least 2"

    # Consider uniqueness of structure to check on attachment positions
    attachment_positions = {match[1] for match in ester_matches}
    if len(attachment_positions) < 2:
        return False, "Acyl chains not placed correctly at 1,2-glycerol positions"

    # Ensure the length criteria of acyl chains by counting carbon nodes on chains
    acyl_chain_count = sum(
        1 for atom in mol.GetAtoms()
        if atom.GetAtomicNum() == 6 and any(bond.GetBondTypeAsDouble() == 1.0 for bond in atom.GetBonds())
    )
    
    if acyl_chain_count < 20:
        return False, "Acyl chains too short to classify"

    return True, "Contains CDP group, appropriate glycerol backbone, and properly attached acyl chains"