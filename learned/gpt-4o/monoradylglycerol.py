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
    # Glycerol backbone pattern (allowing for one substitution): (two OH groups on adjoining carbons)
    glycerol_pattern = Chem.MolFromSmarts("[OX2H][CX4H](O)[CX4](O)O")
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)
    
    if not glycerol_matches:
        return False, "No glycerol backbone found"
    
    # Patterns for acyl (ester linkage), alkyl, and alk-1-enyl groups
    acyl_pattern = Chem.MolFromSmarts("C(=O)[OX2]")
    alkyl_pattern = Chem.MolFromSmarts("C-C")
    alk1_enyl_pattern = Chem.MolFromSmarts("C=C-C")
    
    # Check for exactly one type of lipid chain substituent
    acyl_matches = mol.HasSubstructMatch(acyl_pattern)
    alkyl_matches = mol.HasSubstructMatch(alkyl_pattern)
    alk1_enyl_matches = mol.HasSubstructMatch(alk1_enyl_pattern)
    
    # Count the occurrences
    num_matches = sum([acyl_matches, alkyl_matches, alk1_enyl_matches])
    
    if num_matches != 1:
        return False, f"Expected one substituent, found {num_matches}"
    
    return True, "Contains a glycerol backbone with exactly one lipid substituent"