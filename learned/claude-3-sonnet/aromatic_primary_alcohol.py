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

    # Kekulize the molecule to ensure proper aromatic detection
    Chem.Kekulize(mol, clearAromaticFlags=True)
    # Set aromaticity based on RDKit's model
    Chem.SetAromaticity(mol)

    # Pattern for primary alcohol (-CH2OH) connected to aromatic carbon
    # [cR] matches aromatic carbon in ring
    # [CH2X4] matches CH2 group with 4 connections
    # [OX2H1] matches OH group
    patterns = [
        Chem.MolFromSmarts("[cR][CH2X4][OX2H1]"),  # Aromatic carbon
        Chem.MolFromSmarts("[nR][CH2X4][OX2H1]"),  # Aromatic nitrogen
    ]
    
    matches = []
    for pattern in patterns:
        if pattern is not None:
            matches.extend(mol.GetSubstructMatches(pattern))
    
    if not matches:
        return False, "No primary alcohol (-CH2OH) attached to aromatic ring found"
    
    # Verify at least one valid match
    valid_match = False
    for match in matches:
        aromatic_atom = mol.GetAtomWithIdx(match[0])
        ch2_atom = mol.GetAtomWithIdx(match[1])
        oh_atom = mol.GetAtomWithIdx(match[2])
        
        # Additional checks for true aromaticity
        if not aromatic_atom.GetIsAromatic():
            continue
            
        # Check if the aromatic atom is part of a valid aromatic ring
        ring_info = mol.GetRingInfo()
        if not ring_info.IsAtomInRingOfSize(aromatic_atom.GetIdx(), 5) and \
           not ring_info.IsAtomInRingOfSize(aromatic_atom.GetIdx(), 6):
            continue
            
        # Verify the ring is truly aromatic (all atoms in the ring are aromatic)
        rings = ring_info.AtomRings()
        for ring in rings:
            if aromatic_atom.GetIdx() in ring:
                if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                    # Verify CH2 group
                    if ch2_atom.GetTotalNumHs() == 2:
                        # Verify OH group
                        if oh_atom.GetTotalNumHs() >= 1:
                            valid_match = True
                            break
        
        if valid_match:
            break
    
    if not valid_match:
        return False, "No valid aromatic primary alcohol structure found"

    return True, "Contains primary alcohol (-CH2OH) attached to aromatic ring"