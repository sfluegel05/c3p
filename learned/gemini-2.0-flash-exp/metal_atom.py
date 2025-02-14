"""
Classifies: CHEBI:33521 metal atom
"""
from rdkit import Chem
from rdkit.Chem import PeriodicTable

def is_metal_atom(smiles: str):
    """
    Determines if a molecule is a metal atom based on its SMILES string.
    A metal atom is defined as a single atom of a metallic element.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        tuple(bool, str): (True, reason) if molecule is a metal atom, (False, reason) otherwise.
                          If parsing fails return (None, None)
    """
    try:
        # Parse SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False, "Invalid SMILES string"

        # Check if the molecule is a single atom
        if mol.GetNumAtoms() != 1:
            return False, "Not a single atom."

        atom = mol.GetAtomWithIdx(0)
        symbol = atom.GetSymbol()

        # Remove isotopic number (if present) from the symbol to extract the element's symbol
        element_symbol = ''.join(filter(str.isalpha, symbol))


        pt = PeriodicTable.GetPeriodicTable()
        atomic_num = pt.GetAtomicNumber(element_symbol)
        if atomic_num == 0:
          return False, f"Invalid atom: {element_symbol}"

        element_symbol = pt.GetElementSymbol(atomic_num)
       

        # Check if the atom is a metal
        if atom.IsMetal():
             return True, "Single metal atom"
        else:
            return False, f"Single atom, but not a metal: {element_symbol}"


    except Exception as e:
        return None, None