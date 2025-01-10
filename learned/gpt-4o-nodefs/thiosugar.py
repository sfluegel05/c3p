"""
Classifies: CHEBI:73754 thiosugar
"""
from rdkit import Chem

def is_thiosugar(smiles: str):
    """
    Determines if a molecule is a thiosugar based on its SMILES string.
    A thiosugar is characterized by a sugar-like cyclic structure where
    one of the oxygen atoms is replaced by a sulfur atom.

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

    # Define more comprehensive thiosugar patterns
    # Using patterns more representative of the actual thiosugars known
    thiosugar_patterns = [
        Chem.MolFromSmarts("[C@H]1(O[C@H]([C@@H](O)[C@H](O)[C@H]1O)CO)S"),  # Example: Replacing an oxygen in pyranose rings
        Chem.MolFromSmarts("OC(C)SC1OC(CO)C(O)C(O)C1O")  # A specific thioglycoside pattern 
    ]

    # Check for thiosugar pattern match
    has_thiosugar = any(mol.HasSubstructMatch(pattern) for pattern in thiosugar_patterns)
    if not has_thiosugar:
        return False, "No thiosugar pattern found in the molecule"
    
    # Further confirm by checking essential sugar-like properties
    num_oh_groups = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'O' and any(n.GetSymbol() == 'H' for n in atom.GetNeighbors()))
    
    if num_oh_groups < 3:
        return False, "Insufficient OH groups for a sugar-like structure, found only {}".format(num_oh_groups)

    return True, "Molecule is a thiosugar: Satisfies pattern and functional group requirements"

__metadata__ = {
    'chemical_class': {   
        'name': 'thiosugar',
        'definition': 'Sugars in which one of the oxygen atoms is replaced by sulfur'
    }
}