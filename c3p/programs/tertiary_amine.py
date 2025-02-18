"""
Classifies: CHEBI:32876 tertiary amine
"""
"""
Classifies: tertiary amine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_tertiary_amine(smiles: str):
    """
    Determines if a molecule is a tertiary amine based on its SMILES string.
    A tertiary amine is a nitrogen atom bonded to three carbon atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tertiary amine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get all N atoms
    nitrogen_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]

    if not nitrogen_atoms:
      return False, "No nitrogen atoms found"
    
    # check if each N is bonded to exactly three C's
    for atom in nitrogen_atoms:
        carbon_neighbors = 0
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:
                carbon_neighbors += 1
        
        if carbon_neighbors == 3:
          return True, "Tertiary amine"
    
    return False, "Not a tertiary amine"