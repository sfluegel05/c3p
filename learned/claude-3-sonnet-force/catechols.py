"""
Classifies: CHEBI:33566 catechols
"""
"""
Classifies: CHEBI:38181 catechol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_catechols(smiles: str):
    """
    Determines if a molecule contains a catechol component based on its SMILES string.
    A catechol component is defined as an o-diphenol group, potentially with additional substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule contains a catechol component, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a broader catechol pattern
    catechol_pattern = Chem.MolFromSmarts("[c;r5]1:c(:c(:c(:c(:c:1)-[OH]):c-[OH])-*)-*")
    
    # Check if the molecule contains the catechol pattern
    catechol_matches = mol.GetSubstructMatches(catechol_pattern)
    if catechol_matches:
        # Additional checks or heuristics can be added here
        # to ensure the identified substructure is a valid catechol component
        return True, "Contains a catechol component"
    else:
        return False, "Does not contain a catechol component"