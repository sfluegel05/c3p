"""
Classifies: CHEBI:33857 aromatic primary alcohol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_aromatic_primary_alcohol(smiles: str):
    """
    Determines if a molecule is an aromatic primary alcohol based on its SMILES string.
    An aromatic primary alcohol has a -CH2OH group directly attached to an aromatic ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aromatic primary alcohol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for any aromatic atoms
    if not any(atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return False, "No aromatic rings found"

    # Pattern for primary alcohol (-CH2OH) connected to aromatic carbon
    # [c,n,o,s] matches any aromatic atom
    # [CH2X4] matches CH2 group with 4 connections
    # [OX2H1] matches OH group
    pattern = Chem.MolFromSmarts("[c,n,o,s][CH2X4][OX2H1]")
    
    if not mol.HasSubstructMatch(pattern):
        # Try alternate pattern for cases where H is not explicit
        pattern2 = Chem.MolFromSmarts("[c,n,o,s][CH2X4][OX2]")
        if not mol.HasSubstructMatch(pattern2):
            return False, "No primary alcohol (-CH2OH) attached to aromatic ring found"
    
    # Get matches to analyze the environment
    matches = mol.GetSubstructMatches(pattern) or mol.GetSubstructMatches(pattern2)
    
    # Verify at least one valid match
    valid_match = False
    for match in matches:
        aromatic_atom = mol.GetAtomWithIdx(match[0])
        ch2_atom = mol.GetAtomWithIdx(match[1])
        oh_atom = mol.GetAtomWithIdx(match[2])
        
        # Verify aromatic atom is part of an aromatic ring
        if aromatic_atom.GetIsAromatic():
            # Verify CH2 group
            if ch2_atom.GetTotalNumHs() == 2:
                # Verify OH group
                if oh_atom.GetTotalNumHs() >= 1:
                    valid_match = True
                    break
    
    if not valid_match:
        return False, "No valid aromatic primary alcohol structure found"

    return True, "Contains primary alcohol (-CH2OH) attached to aromatic ring"