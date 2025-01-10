"""
Classifies: CHEBI:76575 monoradylglycerol
"""
from rdkit import Chem

def is_monoradylglycerol(smiles: str):
    """
    Determines if a molecule is a monoradylglycerol based on its SMILES string.
    A monoradylglycerol is a derivative of glycerol where one hydroxyl group is esterified 
    with a fatty acid or similar aliphatic or aromatic acid.

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

    # Broader SMARTS pattern for mostly intact glycerol with one esterified group
    glycerol_ester_pattern = Chem.MolFromSmarts("O[C@H]([CX4H2,CX3H1]!@)O[C](=O)C")  # Allow flexibility in ester linkage
    if not mol.HasSubstructMatch(glycerol_ester_pattern):
        return False, "No suitable esterified glycerol backbone found"

    # Validate presence of just a single ester linkage on glycerol
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester linkages, expected exactly 1"

    # Extended recognition of carbon chains attached to the ester
    varied_hydrocarbon_chain_pattern = Chem.MolFromSmarts("C[CH2,C](C)[C,C]=C")
    hydrocarbon_matches = mol.GetSubstructMatches(varied_hydrocarbon_chain_pattern)
    if len(hydrocarbon_matches) < 1:
        return False, "Missing or too short hydrocarbon chain"

    # Check for any specific stereochemistry for sn-configurations
    sn_glycerol_pattern = Chem.MolFromSmarts("C[C@@H](CO)OC")
    stereo_info = ""
    if mol.HasSubstructMatch(sn_glycerol_pattern):
        stereo_info = "sn-glycerol configuration detected"
    else:
        stereo_info = "Standard glycerol detected"

    return True, f"Contains a glycerol backbone with one esterified hydroxyl group and a hydrocarbon chain. {stereo_info}"