"""
Classifies: CHEBI:71971 neoflavonoid
"""
from rdkit import Chem

def is_neoflavonoid(smiles: str):
    """
    Determines if a molecule is a neoflavonoid based on its SMILES string.
    A neoflavonoid is any 1-benzopyran with an aryl substituent at position 4.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a neoflavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES 
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define enhanced benzopyran pattern and position-4 aryl attachment
    benzopyran_pattern = Chem.MolFromSmarts("c1ccc2occc2c1")  # Refined for benzopyran
    aryl_pattern = Chem.MolFromSmarts("c1ccc[cH]c1")  # Matching phenyl/aryl groups

    # Find the benzopyran in the molecule
    benzopyran_matches = mol.GetSubstructMatches(benzopyran_pattern)
    if not benzopyran_matches:
        return False, "No 1-benzopyran structure found"

    # Check each benzopyran structure for an aryl group at position 4
    for match in benzopyran_matches:
        # Assume the fourth atom in the match is indeed position 4 (considering connectivity)
        pos_4_atom_idx = match[3]
        atom = mol.GetAtomWithIdx(pos_4_atom_idx)
        
        # Checking if any neighbor of pos_4_atom_idx is an aryl group
        for neighbor in atom.GetNeighbors():
            # Check for aryl group attachment
            if mol.HasSubstructMatch(aryl_pattern, atoms=[neighbor.GetIdx()]):
                return True, "Aryl substituent found at position 4"
    
    return False, "No aryl substitution at position 4 in the 1-benzopyran"