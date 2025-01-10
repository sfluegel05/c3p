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

    # Pattern for primary alcohol (-CH2OH) connected to aromatic atom
    # [a] matches any aromatic atom
    # [CH2X4] matches CH2 group with exactly 4 connections
    # [OH1] matches hydroxyl group with exactly 1 hydrogen
    pattern = Chem.MolFromSmarts("[a][CH2X4][OH1]")
    
    if not pattern or not mol.HasSubstructMatch(pattern):
        return False, "No primary alcohol (-CH2OH) attached to aromatic ring found"
    
    # Get all matches
    matches = mol.GetSubstructMatches(pattern)
    
    # Verify each match
    for match in matches:
        aromatic_atom = mol.GetAtomWithIdx(match[0])
        ch2_atom = mol.GetAtomWithIdx(match[1])
        oh_atom = mol.GetAtomWithIdx(match[2])
        
        # Verify aromaticity
        if not aromatic_atom.GetIsAromatic():
            continue
            
        # Get the ring that contains the aromatic atom
        rings = mol.GetRingInfo().AtomRings()
        for ring in rings:
            if aromatic_atom.GetIdx() in ring:
                # Check if the ring is aromatic (all atoms in ring are aromatic)
                ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
                if all(atom.GetIsAromatic() for atom in ring_atoms):
                    # Verify CH2OH group
                    if (ch2_atom.GetTotalNumHs() == 2 and 
                        oh_atom.GetTotalNumHs() == 1 and
                        ch2_atom.GetHybridization() == Chem.HybridizationType.SP3):
                        return True, "Contains primary alcohol (-CH2OH) attached to aromatic ring"
    
    return False, "No valid aromatic primary alcohol structure found"