"""
Classifies: CHEBI:76578 diradylglycerol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_diradylglycerol(smiles: str):
    """
    Determines if a molecule is a diradylglycerol based on its SMILES string.
    A diradylglycerol is a glycerol backbone with two fatty acid chains attached via ester bonds or ethers.

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
    
    # Define a glycerol backbone pattern: Carbon chain COC or COCCO 
    glycerol_pattern = Chem.MolFromSmarts("C(CO)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Identify two ester or ether groups
    ester_or_ether_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H0] | [OX2H0]")
    ester_ether_matches = mol.GetSubstructMatches(ester_or_ether_pattern)
    if len(ester_ether_matches) < 2:
        return False, f"Found {len(ester_ether_matches)} ester or ether groups, need at least 2"
    
    return True, "Contains a glycerol backbone with two fatty acid chains attached via ester bonds or ethers"