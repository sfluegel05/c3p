"""
Classifies: CHEBI:73754 thiosugar
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_thiosugar(smiles: str):
    """
    Determines if a molecule is a thiosugar based on its SMILES string.
    A thiosugar is characterized by a sugar-like cyclic structure where
    one of the oxygen atoms has been replaced by a sulfur atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a thiosugar, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the thiosugar pattern as a cyclic structure with sulfur
    thiosugar_pattern = Chem.MolFromSmarts("C1CS[C@H]1")  # Simple representation of a ring containing carbon and sulfur

    if not mol.HasSubstructMatch(thiosugar_pattern):
        return False, "No thiosugar pattern (ring with sulfur) found"
    
    # Check for a minimum number of OH groups to confirm sugar-like properties
    num_oh_groups = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'O' and any(n.GetSymbol() == 'H' for n in atom.GetNeighbors()))
    
    if num_oh_groups < 3:
        return False, f"Insufficient OH groups for sugar-like properties, found {num_oh_groups}"
    
    return True, "Contains thiosugar pattern with sulfur in sugar-like cyclic structure"


__metadata__ = {   
    'chemical_class': {   
        'name': 'thiosugar',
        'definition': 'Sugars in which one of the oxygen atoms is replaced by sulfur'
    }
}