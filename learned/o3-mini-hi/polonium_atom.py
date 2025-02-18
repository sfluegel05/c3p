"""
Classifies: CHEBI:33313 polonium atom
"""
"""
Classifies: Polonium atom (radioactive metallic element, discovered 1898 by Marie Curie, named after Poland)
Examples: [218Po], [202Po], etc.
"""
from rdkit import Chem

def is_polonium_atom(smiles: str):
    """
    Determines if a molecule is a polonium atom of an allowed isotope based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is a valid polonium atom with an allowed isotope, False otherwise.
        str: Reason for classification.
    """
    
    # Parse the SMILES string using rdkit
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # For a single atom, the molecule should have exactly one atom.
    if mol.GetNumAtoms() != 1:
        return False, "The SMILES does not represent a single-atom molecule"
    
    atom = mol.GetAtomWithIdx(0)
    
    # Check if the atom is polonium (atomic number 84)
    if atom.GetAtomicNum() != 84:
        return False, "The atom is not polonium"
    
    # Get the isotope number from the atom (0 if not specified)
    isotope = atom.GetIsotope()
    if isotope == 0:
        return False, "No isotope specified for the polonium atom"
    
    # Allowed isotopes based on provided examples
    allowed_isotopes = {218, 202, 208, 209, 190, 201, 216, 203, 210, 
                        217, 206, 211, 193, 194, 195, 197, 204, 200, 
                        196, 199, 191, 198, 214}
    
    if isotope not in allowed_isotopes:
        return False, f"Isotope {isotope} is not recognized as a valid polonium isotope in this classification"
    
    return True, f"Valid polonium atom with isotope {isotope}"