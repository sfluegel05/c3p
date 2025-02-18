"""
Classifies: CHEBI:33313 polonium atom
"""
"""
Classifies: Polonium Atom
A radioactive metallic element discovered in 1898 by Marie Sklodowska Curie, 
named after her home country Poland (Latin Polonia).
This function checks if the input SMILES string represents a polonium atom 
with an explicit isotopic mass specifier.
"""
from rdkit import Chem

def is_polonium_atom(smiles: str):
    """
    Determines if a molecule is a polonium atom with an isotope mass number based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule (e.g., "[210Po]")

    Returns:
        bool: True if the SMILES represents a polonium atom with an isotope, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Molecule must consist of exactly one atom
    atoms = mol.GetAtoms()
    if len(atoms) != 1:
        return False, f"Expected a single atom but found {len(atoms)} atoms"
    
    atom = atoms[0]
    
    # Check if the atom symbol is Polonium (Po)
    if atom.GetSymbol() != "Po":
        return False, f"Atom is {atom.GetSymbol()}, not Polonium (Po)"
    
    # Check that an isotope (mass number) is provided, e.g., 190 in [190Po].
    # The GetIsotope() method returns 0 when no isotope is specified.
    isotope = atom.GetIsotope()
    if not isotope:
        return False, "No isotope specified; expected a specific polonium isotope (e.g., [210Po])"
    
    # Optionally, you could verify the isotope number is within a reasonable range.
    # For now, we simply confirm that a mass number is given.
    return True, f"Valid polonium atom with isotope mass: {isotope}"