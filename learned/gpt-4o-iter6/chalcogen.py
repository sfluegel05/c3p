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
    if mol.GetNumAtoms() == 1:
        # Get the single atom in the molecule
        atom = mol.GetAtomWithIdx(0)
        
        # Define atomic numbers for chalcogen elements
        chalcogen_atomic_numbers = {8, 16, 34, 52, 84}  # O, S, Se, Te, Po
        
        # Check if the atom is a chalcogen
        if atom.GetAtomicNum() in chalcogen_atomic_numbers and atom.GetFormalCharge() == 0:
            return True, f"Contains neutral chalcogen element with atomic number: {atom.GetAtomicNum()}"
        else:
            return False, "No neutral chalcogen elements found"
    else:
        return False, "SMILES does not represent a single atom or represents a molecule"