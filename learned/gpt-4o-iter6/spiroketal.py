"""
Classifies: CHEBI:72600 spiroketal
"""
from rdkit import Chem

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

    # Find all ring systems
    ring_info = mol.GetRingInfo()
    
    # Iterate through all atoms to find spiro atoms
    spiro_atoms = []
    for atom in mol.GetAtoms():
        if ring_info.NumAtomRings(atom.GetIdx()) == 2:
            # Atom is part of exactly 2 rings
            spiro_atoms.append(atom)
    
    if not spiro_atoms:
        return False, "No spiro atom found"

    # Define a pattern that matches a ketal structure 
    ketal_pattern = Chem.MolFromSmarts("[C]([O])([O])")
    # Go through each spiro atom and verify if it forms the ketal linkage
    for spiro_atom in spiro_atoms:
        if mol.HasSubstructMatch(ketal_pattern, recalculate=True):
            return True, "Found a spiroketal with spiro atom forming ketal linkage"
    
    return False, "No ketal linkage found on spiro atom"

# Example usage
smiles_examples = [
    "O=C1O[C@@H]2C[C@H](O[C@@H](C2)CC=C(C[C@H](C=CC=C([C@]4[C@H]1C=C([C@@H](O)C4O)C)O)C)C)C",
    "CC[C@H]([C@H]1CC[C@H](C)[C@@H](O1)[C@@H](C)[C@H](O)[C@H](C)C(=O)[C@H](CC)[C@H]1O",
]
for smiles in smiles_examples:
    result, reason = is_spiroketal(smiles)
    print(f"SMILES: {smiles} - Is Spiroketal: {result}, Reason: {reason}")