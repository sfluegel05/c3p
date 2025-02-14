"""
Classifies: CHEBI:1722 3beta-hydroxy-Delta(5)-steroid
"""
from rdkit import Chem

def is_3beta_hydroxy_Delta_5__steroid(smiles: str):
    """
    Determines if a molecule is a 3beta-hydroxy-Delta(5)-steroid based on its SMILES string.
    A 3beta-hydroxy-Delta(5)-steroid is a 3beta-hydroxy-steroid that contains a double bond
    between positions 5 and 6.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 3beta-hydroxy-Delta(5)-steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Simplified 3beta-hydroxy group pattern
    hydroxy_pattern = Chem.MolFromSmarts("[C@@H]([O])[C]")

    # Steroid nucleus pattern (ABCD rings)
    steroid_pattern = Chem.MolFromSmarts("[#6]1-[#6]2-[#6]3-[#6]4-[#6]([#6][#6]3)-[#6]([#6]2)-[#6]1") 

    # Delta(5) double bond pattern to match any double bond in rings
    delta5_pattern = Chem.MolFromSmarts("[C]=[C]")

    # Check for 3beta-hydroxy group
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "No 3beta-hydroxy group match found."

    # Check for steroid backbone
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone detected."

    # Check for Delta(5) double bond within a ring
    delta5_matches = mol.GetSubstructMatches(delta5_pattern)
    delta5_within_ring = any(mol.GetRingInfo().NumAtomRings(match[0]) > 0 and mol.GetRingInfo().NumAtomRings(match[1]) > 0 for match in delta5_matches)
    
    if not delta5_within_ring:
        return False, "No specific Delta(5) double bond match found in ring."

    return True, "Molecule is classified as 3beta-hydroxy-Delta(5)-steroid."