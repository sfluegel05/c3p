"""
Classifies: CHEBI:137419 secondary ammonium ion
"""
"""
Classifies: CHEBI:35275 secondary ammonium ion
"""
from rdkit import Chem

def is_secondary_ammonium_ion(smiles: str):
    """
    Determines if a molecule is a secondary ammonium ion based on its SMILES string.
    A secondary ammonium ion is obtained by protonation of a secondary amine,
    resulting in a nitrogen atom with a +1 formal charge bonded to two substituents
    (any atoms) and two hydrogen atoms.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a secondary ammonium ion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Initialize a flag to indicate presence of secondary ammonium ion
    has_secondary_ammonium = False
    
    # Iterate over all atoms in the molecule
    for atom in mol.GetAtoms():
        # Check if the atom is nitrogen with a +1 formal charge
        if atom.GetAtomicNum() == 7 and atom.GetFormalCharge() == 1:
            # Get the total number of hydrogens attached to the nitrogen (implicit and explicit)
            num_H = atom.GetTotalNumHs()
            # Get the degree (number of neighboring atoms)
            degree = atom.GetDegree()
            # Check that nitrogen has two hydrogens and is connected to two other atoms
            if num_H == 2 and degree == 4:
                # Count the number of non-hydrogen neighbors
                non_H_neighbors = sum(1 for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() > 1)
                # Check if there are exactly two non-hydrogen neighbors
                if non_H_neighbors == 2:
                    has_secondary_ammonium = True
                    break  # Stop after finding one instance
    
    if has_secondary_ammonium:
        return True, "Contains a secondary ammonium ion group (protonated secondary amine)"
    else:
        return False, "Does not contain a secondary ammonium ion group"