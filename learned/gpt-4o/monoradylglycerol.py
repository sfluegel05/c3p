"""
Classifies: CHEBI:76575 monoradylglycerol
"""
from rdkit import Chem

def is_monoradylglycerol(smiles: str):
    """
    Determines if a molecule is a monoradylglycerol based on its SMILES string.
    A monoradylglycerol has a glycerol backbone with a single acyl, alkyl, or alk-1-enyl 
    substituent at an unspecified position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoradylglycerol, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string 
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Identify glycerol backbone with one substituent
    # Glycerol backbone pattern: (two OH groups on adjoining carbons)
    glycerol_pattern = Chem.MolFromSmarts("C(CO)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone with two OH groups found"
    
    # Patterns for acyl, alkyl, and alk-1-enyl linked to glycerol
    acyl_pattern = Chem.MolFromSmarts("[$([CX3](=O)[OX2H1]),$([CX3](=O)[O][CX4])]")
    ether_pattern = Chem.MolFromSmarts("[OX2][CX4]")
    vinyl_ether_pattern = Chem.MolFromSmarts("[C,C]=[C][O]")
    
    # Count the occurrence of these patterns
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)
    ether_matches = mol.GetSubstructMatches(ether_pattern)
    vinyl_ether_matches = mol.GetSubstructMatches(vinyl_ether_pattern)
    
    # Check if it contains exactly one type of lipid chain substituent
    num_matches = len(acyl_matches) + len(ether_matches) + len(vinyl_ether_matches)
    
    if num_matches != 1:
        return False, f"Expected one lipid substituent, found {num_matches}"
    
    return True, "Contains a glycerol backbone with one lipid substituent"