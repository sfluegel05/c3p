"""
Classifies: CHEBI:72600 spiroketal
"""
from rdkit import Chem

def is_spiroketal(smiles: str):
    """
    Determines if a molecule is a spiroketal based on its SMILES string.
    A spiroketal is defined by a ketal carbon that is the only common atom of two rings.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a spiroketal, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES string"

    # SMARTS pattern for a spiroketal: A carbon single-bonded to two oxygens, each part of a distinct ring
    spiroketal_pattern = Chem.MolFromSmarts("[C](O)(O)@[R1]@[R2]")  # Simplified for demonstration

    # Check for spiroketal substructure in the molecule
    if not spiroketal_pattern or not mol.HasSubstructMatch(spiroketal_pattern):
        return False, "No spiroketal structure found"

    # Verify that this setup forms a valid spiro configuration with two distinct rings
    matches = mol.GetSubstructMatches(spiroketal_pattern)
    for match in matches:
        ketal_carbon = match[0]
        rings = [set(ring) for ring in Chem.GetSymmSSSR(mol) if ketal_carbon in ring]
        if len(set.union(*rings)) != len(rings[0].union(rings[1])) or len(rings) != 2:
            return False, "No valid spiroketal: ketal carbon shared exactly by two rings"
    
    return True, "Contains ketal carbon common to two rings (spiroketal structure)"