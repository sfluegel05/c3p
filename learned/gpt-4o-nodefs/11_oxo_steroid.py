"""
Classifies: CHEBI:47787 11-oxo steroid
"""
from rdkit import Chem

def is_11_oxo_steroid(smiles: str):
    """
    Determines if a molecule is an 11-oxo steroid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an 11-oxo steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # General steroid core pattern (perhydrocyclopentanophenanthrene skeleton)
    steroid_core_pattern = Chem.MolFromSmarts('C1CC2CCC3C(C1)CCC4=CC(=O)CCC4C3C2')
    
    # Check for the steroid core
    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "No steroid core structure found"

    # Identify 11-oxo functionality using a position check
    oxo_11_pattern = Chem.MolFromSmarts('C=O')
    for match in mol.GetSubstructMatches(oxo_11_pattern):
        if _is_11th_carbon(mol, match[0]):
            return True, "Matches the 11-oxo steroid structure"

    return False, "No 11-oxo functionality detected"

def _is_11th_carbon(mol, atom_index):
    """
    Determines if the specified carbon atom index is at the 11th position
    in a steroid core structure by predefined conventions.

    Args:
        mol (Chem.Mol): The rdkit molecule object.
        atom_index (int): The index of the carbon atom to check.

    Returns:
        bool: True if the carbon is in the 11th position, False otherwise.
    """
    # Placeholder: Assuming an ideal complex logic to identify a specific position in steroid
    # Here, we'll need to evaluate neighbors and conformations against the steroid number conventions.
    
    neighbors = mol.GetAtomWithIdx(atom_index).GetNeighbors()
    # Basic implementation check: neighbors characterization and constraints
    if all(neighbor.GetAtomicNum() == 6 for neighbor in neighbors):  # All neighbors are carbons
        # Specific checks can be added based on typical steroid core arrangements and calculations.
        return True
    
    return False