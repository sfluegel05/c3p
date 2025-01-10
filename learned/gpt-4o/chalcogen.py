"""
Classifies: CHEBI:33303 chalcogen
"""
from rdkit import Chem

def is_chalcogen(smiles: str):
    """
    Determines if a SMILES string represents a single chalcogen atom, including isotopes.
    
    Args:
        smiles (str): SMILES string of a chemical entity

    Returns:
        bool: True if the SMILES represents a single chalcogen atom or its isotope, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # There should be exactly one atom in the molecule for it to be a standalone chalcogen element
    if mol.GetNumAtoms() != 1:
        return False, "SMILES represents a compound, not a single atom"

    # Get the single atom in the molecule
    atom = mol.GetAtomWithIdx(0)
    
    # Consider atomic numbers for chalcogens: O (8), S (16), Se (34), Te (52), Po (84)
    chalcogen_atomic_numbers = {8, 16, 34, 52, 84}
    
    # Check if the atom is a chalcogen
    if atom.GetAtomicNum() in chalcogen_atomic_numbers:
        # Allow isotopes, which have the same atomic number
        return True, f"The atom is a chalcogen: {atom.GetSymbol()}"

    return False, "No chalcogen atom present"