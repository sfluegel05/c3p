"""
Classifies: CHEBI:35915 sterol ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_sterol_ester(smiles: str):
    """
    Determines if a molecule is a sterol ester based on its SMILES string.
    A sterol ester is a steroid with an ester linkage formed by condensation of a carboxylic acid
    with the 3-hydroxy group of a sterol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sterol ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Flexibly recognize the steroid skeleton (cyclopentanoperhydrophenanthrene core)
    steroid_pattern = Chem.MolFromSmarts("C1CCC2C(C1)CCC3C2CCC4=C3CCCC4")

    # Recognize ester linkage (C(=O)O) but focus on it being attached in the right way
    ester_pattern = Chem.MolFromSmarts("C(=O)O")

    # Check for steroid skeleton
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid skeleton found"

    # Check for ester linkage
    ester_matches = mol.GetSubstructMatches(ester_pattern)

    if not ester_matches:
        return False, "No ester linkage found"

    # Further refine to ensure 3-hydroxy specific attachment
    for match in ester_matches:
        # Check if the ester is connected to a possible secondary hydroxyl position in steroids
        # This step can be more chemistry aware if needed for more accurate filtering
        ester_atom_index = matches[0]  # Index of O in the ester from pattern match
        connected_atoms = mol.GetAtomWithIdx(ester_atom_index).GetNeighbors()

        # Ensure one connected atom is part of the steroid core
        if any(mol.GetSubstructMatch(steroid_pattern).count(atom.GetIdx()) > 0 for atom in connected_atoms):
            return True, "Contains a steroid skeleton with ester linkage"

    return False, "Ester linkage not found in the appropriate 3-hydroxy position"

# Test the function with a sterol ester example
smiles_example = "CCCCCCCCCCCCCCCCC(=O)O[C@H]1CC[C@]2(C)[C@H]3CC[C@]4(C)[C@H](CC[C@H]4[C@@H]3CC=C2C1)[C@H](C)CCCC(C)C"
print(is_sterol_ester(smiles_example))