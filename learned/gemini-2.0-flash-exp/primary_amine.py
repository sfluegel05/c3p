"""
Classifies: CHEBI:32877 primary amine
"""
from rdkit import Chem

def is_primary_amine(smiles: str):
    """
    Determines if a molecule is a primary amine based on its SMILES string.
    A primary amine is defined as a nitrogen atom bonded to two hydrogen atoms and one other substituent.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a primary amine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES, remove any salt info
    smiles = smiles.split(".")[0]
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check each nitrogen atom
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7: # Check if the atom is Nitrogen
            num_h = 0
            num_non_h = 0
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 1: # If neighbor is H, then count
                    num_h +=1
                else:
                   num_non_h +=1 #if it is not H, then count it as non-H
            if num_non_h == 1 and (num_h <= 2):
                return True, "Contains a primary amine"

    return False, "Does not contain a primary amine"