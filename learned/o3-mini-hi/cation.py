"""
Classifies: CHEBI:36916 cation
"""
"""
Classifies chemical entities of the class cation:
A monoatomic or polyatomic species having one or more elementary charges of the proton.
"""

from rdkit import Chem

def is_cation(smiles: str):
    """
    Determines if a molecule is a cation based on its SMILES string.
    A cation must have a net formal charge greater than zero.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is a cation (net positive charge), False otherwise.
        str: Explanation of the result.
    """
    # Parse the SMILES string to create a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Calculate the net formal charge by summing the formal charges of all atoms
    net_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    
    # Check if the net charge is greater than zero (cation condition)
    if net_charge > 0:
        return True, f"Molecule has a net positive charge of {net_charge}"
    else:
        return False, f"Molecule has a net charge of {net_charge} (not a cation)"