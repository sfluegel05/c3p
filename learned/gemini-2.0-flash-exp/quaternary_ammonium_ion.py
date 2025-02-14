"""
Classifies: CHEBI:35267 quaternary ammonium ion
"""
from rdkit import Chem

def is_quaternary_ammonium_ion(smiles: str):
    """
    Determines if a molecule is a quaternary ammonium ion based on its SMILES string.
    A quaternary ammonium ion is a positively charged nitrogen with four single bonds to non-hydrogen atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a quaternary ammonium ion, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Iterate through atoms to find nitrogens
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7: # Nitrogen has atomic number 7
            # Check if positively charged
            if atom.GetFormalCharge() == 1:
                #Check number of non-hydrogen neighbors
                non_h_neighbors = 0
                for neighbor in atom.GetNeighbors():
                  if neighbor.GetAtomicNum() != 1:
                    non_h_neighbors += 1
                if non_h_neighbors == 4:
                    return True, "Positively charged nitrogen with four bonds to non-hydrogen atoms found."
                
    return False, "No positively charged nitrogen with four bonds to non-hydrogen atoms found."