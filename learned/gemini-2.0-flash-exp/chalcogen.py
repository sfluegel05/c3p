"""
Classifies: CHEBI:33303 chalcogen
"""
from rdkit import Chem

def is_chalcogen(smiles: str):
    """
    Determines if a molecule is a chalcogen based on its SMILES string.
    Chalcogens are elements of group 16: O, S, Se, Te, Po.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a chalcogen, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if molecule contains a single atom
    if mol.GetNumAtoms() != 1:
         return False, "Molecule is not a single atom."

    # Get the atom from the molecule
    atom = mol.GetAtomWithIdx(0)

    # Get atomic symbol
    symbol = atom.GetSymbol()
    
    # Check if it's a chalcogen
    chalcogens = ["O", "S", "Se", "Te", "Po"]
    if symbol in chalcogens:
        return True, "Is a chalcogen atom"
    else:
        return False, "Not a chalcogen atom"