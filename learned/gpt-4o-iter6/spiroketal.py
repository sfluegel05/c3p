"""
Classifies: CHEBI:72600 spiroketal
"""
from rdkit import Chem

def is_spiroketal(smiles: str):
    """
    Determines if a molecule is a spiroketal based on its SMILES string.
    A spiroketal contains a spiro center connected via oxygen atoms to both rings.

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
    
    # Define a pattern that matches the spiro ketal structure: spiro atom
    # bonded to two oxygens each in different rings
    spiro_ketal_pattern = Chem.MolFromSmarts("[C]([O])[O]")
    
    # Go through each atom and verify if it fits the spiroketal pattern
    for atom in mol.GetAtoms():
        # The atom should be in two different rings
        if ring_info.NumAtomRings(atom.GetIdx()) == 2:
            # Verify that the atom matches pattern for spiro ketal
            submol = Chem.PathToSubmol(mol, [atom.GetIdx()])
            if submol.HasSubstructMatch(spiro_ketal_pattern):
                return True, "Found a spiroketal with the correct linkage"

    return False, "No spiroketal found"

# Example usage with the provided SMILES
smiles_examples = [
    "O=C1O[C@@H]2C[C@H](O[C@@H](C2)CC=C(C[C@H](C=CC=C([C@]4[C@H]1C=C([C@@H](O)C4O)C)O)C)C)C",
    "CC[C@H]([C@H]1CC[C@H](C)[C@@H](O1)[C@@H](C)[C@H](O)[C@H](C)C(=O)[C@H](CC)[C@H]1O"
]

# Evaluating the examples
for smiles in smiles_examples:
    result, reason = is_spiroketal(smiles)
    print(f"SMILES: {smiles} - Is Spiroketal: {result}, Reason: {reason}")