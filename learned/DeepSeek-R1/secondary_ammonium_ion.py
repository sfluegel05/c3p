"""
Classifies: CHEBI:137419 secondary ammonium ion
"""
"""
Classifies: CHEBI:95486 secondary ammonium ion
"""
from rdkit import Chem

def is_secondary_ammonium_ion(smiles: str):
    """
    Determines if a molecule is a secondary ammonium ion based on its SMILES string.
    A secondary ammonium ion has a nitrogen atom with +1 charge, two carbon substituents,
    and a total valence of 4 (indicating two bonds to carbons and two hydrogens).

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
    
    # Iterate over all atoms to find potential nitrogen atoms
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'N' and atom.GetFormalCharge() == 1:
            valence = atom.GetTotalValence()
            degree = atom.GetDegree()  # Number of explicit bonds
            
            # Check valence and degree criteria
            if valence == 4 and degree == 2:
                # Check if both neighbors are carbons
                neighbors = atom.GetNeighbors()
                if all(n.GetAtomicNum() == 6 for n in neighbors):
                    return True, "Contains a nitrogen with +1 charge, two carbon substituents, and valence 4"
    
    # If no matching nitrogen found
    return False, "No secondary ammonium ion detected"