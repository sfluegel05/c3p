"""
Classifies: CHEBI:33303 chalcogen
"""
from rdkit import Chem

def is_chalcogen(smiles: str):
    """
    Determines if a SMILES string represents a chalcogen atom.

    Args:
        smiles (str): SMILES string of a chemical entity

    Returns:
        bool: True if the SMILES represents a chalcogen atom, False otherwise
        str: Reason for the classification
    """
    
    # Parse SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check the number of atoms; we are interested only in isolated elements
    if mol.GetNumAtoms() != 1:
        return False, "SMILES should represent a single element"

    # Get the atomic symbol of the lone atom in the molecule
    atom = mol.GetAtomWithIdx(0)
    atom_symbol = atom.GetSymbol()

    # Define group 16 elements (chalcogens)
    chalcogens = {"O", "S", "Se", "Te", "Po"}

    # Check for charge and radical electrons, indicating it should be in natural state
    if atom.GetFormalCharge() != 0 or atom.GetNumRadicalElectrons() > 0:
        return False, "The atom is a chalcogen but not in its neutral natural form"

    # Determine if the atom is a chalcogen
    if atom_symbol in chalcogens:
        return True, f"The atom is a chalcogen: {atom_symbol}"
    else:
        return False, f"The atom is not a chalcogen: {atom_symbol}"