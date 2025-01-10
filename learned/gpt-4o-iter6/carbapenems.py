"""
Classifies: CHEBI:46633 carbapenems
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_carbapenems(smiles: str):
    """
    Determine if a given SMILES string is a carbapenem.

    Args:
        smiles (str): The SMILES string of the chemical entity.
    
    Returns:
        bool: True if SMILES represents a carbapenem, False otherwise.
        str: Explanation of the reasoning for the classification.
    """
    
    # Define the SMARTS pattern for the core carbapenem structure:
    # It's a simplified assumption of a bicyclic structure with sulfur and beta-lactam
    carbapenem_pattern = Chem.MolFromSmarts("C1CNC2=C1SC(=O)N2")
    # Alternatively, if looking for specific positions and combinations, patterns may be adjusted

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the carbapenem pattern in the molecule
    if mol.HasSubstructMatch(carbapenem_pattern):
        return True, "Matches the carbapenem core structure"
    else:
        return False, "Does not match the carbapenem core structure"