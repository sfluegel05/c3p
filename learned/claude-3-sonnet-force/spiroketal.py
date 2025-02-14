"""
Classifies: CHEBI:72600 spiroketal
"""
"""
Classifies: CHEBI:48597 spiroketal
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_spiroketal(smiles: str):
    """
    Determines if a molecule is a spiroketal based on its SMILES string.
    A spiroketal is a cyclic ketal in which the ketal carbon is the only common atom of two rings.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a spiroketal, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for ketal pattern
    ketal_pattern = Chem.MolFromSmarts("[OX2][CX4]([OX2])([OX2])")
    ketal_matches = mol.GetSubstructMatches(ketal_pattern)
    
    for match in ketal_matches:
        ketal_carbon_idx = match[1]  # Index of the ketal carbon
        ketal_carbon = mol.GetAtomWithIdx(ketal_carbon_idx)
        
        # Get the rings containing the ketal carbon
        ring_info = mol.GetRingInfo()
        ketal_carbon_rings = ring_info.AtomRings()[ketal_carbon_idx]
        
        # Check if the ketal carbon is part of exactly two rings
        if len(ketal_carbon_rings) == 2:
            ring1, ring2 = [set(ring) for ring in ketal_carbon_rings]
            
            # Check if the two rings share only the ketal carbon
            if len(ring1.intersection(ring2)) == 1:
                
                # Check if the rings are of appropriate size (4-8 atoms)
                ring_sizes = [len(ring) for ring in ketal_carbon_rings]
                if all(4 <= size <= 8 for size in ring_sizes):
                    
                    # Check if the ketal oxygen atoms are part of the rings
                    ketal_oxygens = [mol.GetAtomWithIdx(idx) for idx in match if idx != ketal_carbon_idx]
                    if all(atom.IsInRing() for atom in ketal_oxygens):
                        return True, "Contains a spiroketal moiety with two rings sharing a ketal carbon"
    
    return False, "No spiroketal moiety found"