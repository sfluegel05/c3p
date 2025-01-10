"""
Classifies: CHEBI:76575 monoradylglycerol
"""
from rdkit import Chem

def is_monoradylglycerol(smiles: str):
    """
    Determines if a molecule is a monoradylglycerol based on its SMILES string.
    A monoradylglycerol is a derivative of glycerol where one hydroxyl group is esterified 
    with a fatty acid or similar aliphatic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoradylglycerol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for glycerol backbone with a single ester group
    glycerol_ester_pattern = Chem.MolFromSmarts("OCC(O)CO[C](=O)")

    if not mol.HasSubstructMatch(glycerol_ester_pattern):
        return False, "No glycerol backbone with esterification found"

    # Count the ester linkages (should be exactly 1)
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester linkages, need exactly 1"

    # Ensure there is a sufficiently long hydrocarbon chain attached
    # Assuming a minimum chain length of 4 carbons for fatty acids
    hydrocarbon_chain_pattern = Chem.MolFromSmarts("CCCC")
    hydrocarbon_matches = mol.GetSubstructMatches(hydrocarbon_chain_pattern)
    if len(hydrocarbon_matches) < 1:
        return False, "Missing or too short hydrocarbon chain"

    # Check for any other specific stereochemistry or functional groups
    sn_glycerol_pattern = Chem.MolFromSmarts("[C@H](CO[C](=O)")  # Example for sn-glycerol
    if mol.HasSubstructMatch(sn_glycerol_pattern):
        stereo_info = "sn-glycerol configuration detected"
    else:
        stereo_info = "Standard glycerol detected"

    return True, f"Contains a glycerol backbone with one esterified hydroxyl group and a hydrocarbon chain. {stereo_info}"