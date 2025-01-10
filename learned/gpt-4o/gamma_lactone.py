"""
Classifies: CHEBI:37581 gamma-lactone
"""
from rdkit import Chem

def is_gamma_lactone(smiles: str):
    """
    Determines if a molecule is a gamma-lactone based on its SMILES string.
    A gamma-lactone is defined by a five-membered lactone ring, possibly showing various substitution patterns.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a gamma-lactone, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string to create a molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a broader range of SMARTS patterns for 5-membered lactone rings including variations.
    gamma_lactone_patterns = [
        Chem.MolFromSmarts("C1OC(=O)CC1"),      # Simple C=O, basic structure
        Chem.MolFromSmarts("C1C(=O)OC=C1"),     # Inclusion of carbonyl variation
        Chem.MolFromSmarts("C1OC(=O)C=C1"),     # Consider unsaturated variants
        Chem.MolFromSmarts("C1C(=O)COC1")       # Carbonyl in ring variation
    ]
    
    # Check for the presence of any gamma-lactone ring using the defined patterns.
    for pattern in gamma_lactone_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains gamma-lactone ring"
    
    return False, "No gamma-lactone ring found"

__metadata__ = {
   'chemical_class': {
       'name': 'gamma-lactone',
       'definition': 'A lactone having a five-membered lactone ring.',
   }
}