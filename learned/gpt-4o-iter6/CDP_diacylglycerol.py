"""
Classifies: CHEBI:17962 CDP-diacylglycerol
"""
from rdkit import Chem

def is_CDP_diacylglycerol(smiles: str):
    """
    Determines if a molecule is a CDP-diacylglycerol based on its SMILES string.
    A CDP-diacylglycerol has a glycerol backbone with two acyl chains at positions 1 and 2, and a cytidine diphosphate group.

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

    # Look for CDP group: more general approach covering variations
    cdp_pattern = Chem.MolFromSmarts("n1c([C@H]([O])[C@@H]1O)[C@H](O)COPO(=O)(OP(=O)(O)O)O")
    if not mol.HasSubstructMatch(cdp_pattern):
        return False, "No CDP (Cytidine diphosphate) group found"
    
    # Look for glycerol backbone pattern: consider a flexible approach
    glycerol_pattern = Chem.MolFromSmarts("[C@H](CO)[C@@H](O)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found or incorrect stereochemistry match"

    # Look for two ester linkages indicating acyl chains
    ester_pattern = Chem.MolFromSmarts("C(=O)[O][C@H]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} acyl chains, need at least 2"

    # Verify sufficient length of acyl chains by carbon count
    acyl_chain_count = sum(
        1 for atom in mol.GetAtoms()
        if atom.GetAtomicNum() == 6 and len(atom.GetBonds()) == 3 # follows backbone nature
    )

    if acyl_chain_count < 20:
        return False, "Acyl chains too short to classify"

    return True, "Contains proper glycerol, CDP group, and acyl chain lengths consistent with CDP-diacylglycerol"