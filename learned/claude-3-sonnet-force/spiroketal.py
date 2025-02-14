"""
Classifies: CHEBI:72600 spiroketal
"""
"""
Classifies: CHEBI:48597 spiroketal
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

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
    
    # Look for spiro pattern (two rings sharing a single atom)
    spiro_pattern = Chem.MolFromSmarts("[R1,R2]@[R3,R4]")
    spiro_matches = mol.GetSubstructMatches(spiro_pattern)
    if not spiro_matches:
        return False, "No spiro atom found"
    
    # Check if spiro atom is a ketal
    ketal_pattern = Chem.MolFromSmarts("[OX2][CX4]([OX2])([OX2])")
    for match in spiro_matches:
        spiro_atom_idx = match[0]  # Extract the atom index from the match tuple
        spiro_atom = mol.GetAtomWithIdx(spiro_atom_idx)
        if mol.HasSubstructMatch(ketal_pattern, atomIdxList=[spiro_atom.GetIdx()]):
            # Check ring sizes
            ring_info = mol.GetRingInfo()
            ring_systems = ring_info.AtomRings()
            ring_sizes = set()
            for ring in ring_systems:
                if spiro_atom_idx in ring:
                    ring_sizes.add(len(ring))
            if len(ring_sizes) == 2:
                return True, "Contains a spiroketal moiety with two rings sharing a ketal carbon"
    
    return False, "No spiroketal moiety found"