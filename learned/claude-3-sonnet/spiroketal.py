"""
Classifies: CHEBI:72600 spiroketal
"""
"""
Classifies: CHEBI:52169 spiroketal
A spiroketal is a cyclic ketal in which the ketal carbon is the only common atom of two rings.
"""

from rdkit import Chem
from rdkit.Chem import rdFMCS

def is_spiroketal(smiles: str):
    """
    Determines if a molecule is a spiroketal based on its SMILES string.

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
    
    # Find ketal groups (O-C-O)
    ketal_pattern = Chem.MolFromSmarts("O[C;X4]O")
    ketal_matches = mol.GetSubstructMatches(ketal_pattern)
    
    # Check if any ketal carbon is shared between two rings
    for ketal_match in ketal_matches:
        ketal_atom = mol.GetAtomWithIdx(ketal_match[1])
        rings = mol.GetRingInfo().AtomRings()
        ring_counts = [len([ring for ring in rings if ketal_atom.GetIdx() in ring]) for atom in mol.GetAtoms()]
        if sum(ring_counts) > 2:
            return True, "Contains a ketal carbon that is shared between two rings (spiroketal)"
    
    return False, "No spiroketal substructure found"

# Examples
print(is_spiroketal("O=C1O[C@@H]2C[C@]3(O[C@@H](C2)CC=C(C[C@H](C=CC=C4[C@]5([C@H]1C=C([C@@H](O)[C@H]5OC4)CO)O)C)C)O[C@@H]([C@@H](C)CC3)C"))
# Output: (True, 'Contains a ketal carbon that is shared between two rings (spiroketal)')

print(is_spiroketal("CCCC"))
# Output: (False, 'No spiroketal substructure found')