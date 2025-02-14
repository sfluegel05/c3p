"""
Classifies: CHEBI:76578 diradylglycerol
"""
"""
Classifies: CHEBI:24048 diradylglycerol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_diradylglycerol(smiles: str):
    """
    Determines if a molecule is a diradylglycerol based on its SMILES string.
    A diradylglycerol is a glycerol backbone with two substituent groups - either acyl, alkyl, or alk-1-enyl - at any two of the three possible positions.

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

    # Look for glycerol backbone pattern
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Look for ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)

    # Look for ether groups (-O-C-)
    ether_pattern = Chem.MolFromSmarts("[OX2][CX4]")
    ether_matches = mol.GetSubstructMatches(ether_pattern)
    
    # Look for enol ether groups (-O-C=C-)
    enol_ether_pattern = Chem.MolFromSmarts("[OX2][CX3]=[CX3]")
    enol_ether_matches = mol.GetSubstructMatches(enol_ether_pattern)

    total_substituents = len(ester_matches) + len(ether_matches) + len(enol_ether_matches)

    #Check it has exactly two substituents.
    if total_substituents != 2:
        return False, f"Found {total_substituents} substituents, need exactly 2"

    #Check it is not a triglyceride (3 esters)
    if len(ester_matches) == 3:
        return False, "Molecule is a triglyceride, not a diradylglycerol"

    return True, "Contains glycerol backbone with exactly two substituents (acyl, alkyl or alk-1-enyl) attached"