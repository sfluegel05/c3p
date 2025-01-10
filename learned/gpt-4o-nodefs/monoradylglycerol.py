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

    # SMARTS pattern for a glycerol backbone with ester group
    glycerol_ester_pattern = Chem.MolFromSmarts("OCC(O)COC(=O)")

    if not mol.HasSubstructMatch(glycerol_ester_pattern):
        return False, "No glycerol backbone with esterification found"

    # Count the ester linkages (should be exactly 1)
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester linkages, need exactly 1"

    # Match the hydrocarbon chain pattern
    hydrocarbon_chain_pattern = Chem.MolFromSmarts("C[C,C]C")
    hydrocarbon_matches = mol.GetSubstructMatches(hydrocarbon_chain_pattern)
    if len(hydrocarbon_matches) < 1:
        return False, "Missing hydrocarbon chain"

    # Validate the expected glycerol backbone structure
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "Glycerol backbone is not connected properly"

    return True, "Contains a glycerol backbone with one esterified hydroxyl group and a hydrocarbon chain"