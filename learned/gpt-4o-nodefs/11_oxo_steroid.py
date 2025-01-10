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
    
    # General steroid core pattern, potentially more flexible
    # Representing the steroid core as four rings (ABC ring system for steroids)
    steroid_core_pattern = Chem.MolFromSmarts('C1CCC2C3CCC4CCCC(C3)C24C1')
    
    # Check for the steroid core
    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "No steroid core structure found"

    # Identify 11th position ketone (oxo) functionality
    oxo_pattern = Chem.MolFromSmarts('[C]=O')
    ketone_matches = mol.GetSubstructMatches(oxo_pattern)
    
    for match in ketone_matches:
        if _is_11th_carbon_in_steroid_core(mol, match[0]):
            return True, "Matches the 11-oxo steroid structure"

    return False, "No 11-oxo functionality detected"

def _is_11th_carbon_in_steroid_core(mol, atom_index):
    """
    Determines if the specified carbon atom index corresponds to the 11th position
    in typical steroid numbering.

    Args:
        mol (Chem.Mol): The rdkit molecule object.
        atom_index (int): The index of the carbon atom to check.

    Returns:
        bool: True if the carbon is the 11th position, False otherwise.
    """
    # Using more specific checks for connectivity and geometry
    # This assumes the steroid framework has been indexed for reference.
    atom = mol.GetAtomWithIdx(atom_index)
    
    for neighbor in atom.GetNeighbors():
        # Looking for specific connectivity to identify 11th position
        if neighbor.GetAtomicNum() == 6:  # Carbon neighbor
            # Potential further geometric logic to verify specific location in steroid
            ...
    
    # Contextual verification required to ensure the position is correct
    return False  # Implement specific criteria based on steroid schema

# This function implementation is incomplete, and would need specific logic accurately matching steroid number.