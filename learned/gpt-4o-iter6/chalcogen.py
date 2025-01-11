"""
Classifies: CHEBI:33303 chalcogen
"""
from rdkit import Chem

def is_chalcogen(smiles: str):
    """
    Determines if a molecule is a chalcogen element based on its SMILES string.
    Chalcogens are elements belonging to group 16 of the periodic table, including O, S, Se, Te, Po.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a single, neutral chalcogen atom or isotope, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure the molecule consists of a single atom
    if mol.GetNumAtoms() != 1:
        return False, "SMILES does not represent a single atom"

    # Define atomic numbers for chalcogen elements
    chalcogen_atomic_numbers = {8, 16, 34, 52, 84}  # O, S, Se, Te, Po

    # Check if the single atom in the molecule is a chalcogen and is neutral
    atom = mol.GetAtomWithIdx(0)
    atomic_num = atom.GetAtomicNum()
    formal_charge = atom.GetFormalCharge()
    
    # Ensure the atom is a recognized chalcogen and is neutral
    if atomic_num in chalcogen_atomic_numbers and formal_charge == 0:
        return True, f"Contains neutral chalcogen element with atomic number: {atomic_num}"
    
    return False, "No neutral chalcogen elements found"