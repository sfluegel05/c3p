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

    # Define benzopyran pattern with specified position for aryl attachment at position 4
    # Benzopyran with a generic attachment site at position 4 (not specific, just used for matching)
    benzopyran_pattern = Chem.MolFromSmarts("c1cc2ccocc2c1")
    
    # Define an aryl group pattern
    aryl_pattern = Chem.MolFromSmarts("c1cccc(c1)")
    
    # Attempt to match the 1-benzopyran backbone
    benzopyran_matches = mol.GetSubstructMatches(benzopyran_pattern)
    
    if not benzopyran_matches:
        return False, "No 1-benzopyran structure found"

    # Verify position 4 aryl attachment by checking neighbors of atom at position 4
    # In SMARTS indices, assume atom order in pattern as [_, _, _, aryl_position, ...]
    for match in benzopyran_matches:
        pos_4_atom_idx = match[4]  # 4th position in benzopyran pattern
        atom = mol.GetAtomWithIdx(pos_4_atom_idx)
        
        # Check if any neighbor of 4-position atom forms an aryl group
        for neighbor in atom.GetNeighbors():
            if mol.HasSubstructMatch(aryl_pattern, [neighbor.GetIdx()]):
                return True, "Aryl substituent found at position 4"
    
    return False, "No aryl substitution at position 4 in the 1-benzopyran"