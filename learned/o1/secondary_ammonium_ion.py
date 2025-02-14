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
    resulting in a nitrogen atom with a +1 formal charge bonded to two carbon atoms.
    
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
        # Check if the atom is nitrogen with a +1 formal charge and is not aromatic
        if atom.GetAtomicNum() == 7 and atom.GetFormalCharge() == 1 and not atom.GetIsAromatic():
            # Get the neighbors of the nitrogen atom
            neighbors = atom.GetNeighbors()
            # Count the number of carbon neighbors
            carbon_neighbors = sum(1 for neighbor in neighbors if neighbor.GetAtomicNum() == 6)
            # Check if there are exactly two carbon neighbors
            if carbon_neighbors == 2:
                has_secondary_ammonium = True
                break  # Stop after finding one instance
    
    if has_secondary_ammonium:
        return True, "Contains a secondary ammonium ion group (protonated secondary amine)"
    else:
        return False, "Does not contain a secondary ammonium ion group"