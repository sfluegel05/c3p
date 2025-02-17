"""
Classifies: CHEBI:37141 organobromine compound
"""
"""
Classifies: Organobromine compound
Definition: A compound containing at least one carbon–bromine bond.
This version uses a SMARTS matcher to inspect whether any substructure in the molecule 
matches a carbon atom (atomic number 6) directly bonded to a bromine atom.
"""

from rdkit import Chem

def is_organobromine_compound(smiles: str):
    """
    Determines if a molecule is an organobromine compound based on its SMILES string.
    An organobromine compound is defined as a compound containing at least one carbon–bromine bond.
    
    This implementation uses a SMARTS pattern to directly match any bond where a carbon (atomic number 6)
    is connected to a bromine (atomic number 35). This approach is less prone to ambiguity than iterating over
    all bonds for aromatic molecules.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is an organobromine compound, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string into an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS pattern to match a C–Br bond.
    # This pattern will match any bond between a carbon (atomic number 6) and bromine (atomic number 35)
    # regardless of whether the carbon is aromatic or aliphatic.
    pattern = Chem.MolFromSmarts("[#6]-[Br]")
    if pattern is None:
        return False, "Error in SMARTS pattern"
    
    # Test if the molecule has any substructure that matches the C–Br bond pattern.
    if mol.HasSubstructMatch(pattern):
        return True, "Contains at least one carbon–bromine bond"
    else:
        return False, "No carbon–bromine bond found"