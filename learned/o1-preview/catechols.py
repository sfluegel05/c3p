"""
Classifies: CHEBI:33566 catechols
"""
"""
Classifies: catechols
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_catechols(smiles: str):
    """
    Determines if a molecule is a catechol based on its SMILES string.
    A catechol is any compound containing an o-diphenol component, i.e., a benzene ring with two adjacent hydroxyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains a catechol moiety, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define catechol SMARTS pattern (benzene ring with two adjacent hydroxyl groups)
    catechol_pattern = Chem.MolFromSmarts("c1c(O)c(O)cccc1")
    if catechol_pattern is None:
        return False, "Invalid catechol SMARTS pattern"
    
    # Search for catechol pattern in molecule
    if mol.HasSubstructMatch(catechol_pattern):
        return True, "Contains an o-diphenol (catechol) moiety"
    else:
        return False, "No catechol moiety found"
    
__metadata__ = {
   'chemical_class': {
       'id': 'CHEBI:33539',
       'name': 'catechols',
       'definition': 'Any compound containing an o-diphenol component.',
       'parents': []
   }
}