"""
Classifies: CHEBI:76578 diradylglycerol
"""
from rdkit import Chem

def is_diradylglycerol(smiles: str):
    """
    Determines if a molecule is a diradylglycerol based on its SMILES string.
    A diradylglycerol is characterized by a glycerol backbone with two fatty acid chains attached as esters or ethers,
    and one free hydroxyl group (usually at the sn-3 position).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diradylglycerol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Correct glycerol backbone pattern C(CO)CO with potential stereochemistry
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Identify ester: R-C(=O)-O or ether: R-O-R' linkages
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ether_pattern = Chem.MolFromSmarts("COC")
    
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    ether_matches = mol.GetSubstructMatches(ether_pattern)

    # Total should be at least 2 esters/ethers indicating two fatty acid chains
    if len(ester_matches) + len(ether_matches) < 2:
        return False, f"Found {len(ester_matches)} esters and {len(ether_matches)} ethers, need at least two in any combination"
    
    return True, "Contains glycerol backbone with two ester or ether linkages"