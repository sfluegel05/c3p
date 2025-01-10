"""
Classifies: CHEBI:76578 diradylglycerol
"""
from rdkit import Chem

def is_diradylglycerol(smiles: str):
    """
    Determines if a molecule is a diradylglycerol based on its SMILES string.
    A diradylglycerol is a glycerol backbone with two fatty acid chains attached as esters or ethers.

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
    
    # Define glycerol backbone pattern C(CO)CO
    glycerol_pattern = Chem.MolFromSmarts("[C@@H](CO)CO | [C@H](CO)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Identify ester groups: R-C(=O)-O-R' or ether: R-O-R'
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ether_pattern = Chem.MolFromSmarts("COC")
    
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    ether_matches = mol.GetSubstructMatches(ether_pattern)
    
    if len(ester_matches) + len(ether_matches) < 2:
        return False, f"Found {len(ester_matches)} esters and {len(ether_matches)} ethers, need at least two in any combination"
    
    return True, "Contains a glycerol backbone with two ester or ether linkages"