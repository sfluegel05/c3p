"""
Classifies: CHEBI:37143 organofluorine compound
"""
"""
Classifies: CHEBI:35322 organofluorine compound
An organofluorine compound is a compound containing at least one carbon-fluorine bond.
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_organofluorine_compound(smiles: str):
    """
    Determines if a molecule is an organofluorine compound based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organofluorine compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carbon-fluorine bonds using RDKit's built-in functionality
    n_c_f_bonds = Descriptors.Fr_C_F(mol)
    
    # If any carbon-fluorine bond is found, classify as organofluorine compound
    if n_c_f_bonds > 0:
        return True, f"Contains {n_c_f_bonds} carbon-fluorine bond(s)"
    
    # If no carbon-fluorine bond is found, classify as non-organofluorine compound
    return False, "No carbon-fluorine bond found"