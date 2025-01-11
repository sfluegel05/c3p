"""
Classifies: CHEBI:35342 17alpha-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_17alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 17alpha-hydroxy steroid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule is a 17alpha-hydroxy steroid, False otherwise
        str: Reason or structural match explanation
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Determine presence of typical steroid backbone
    steroid_pattern = Chem.MolFromSmarts("C1CCC2C(C1)CCC3C2CCC4(C3(CCC4C4CC(=O)CCC4)C)C")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No typical steroid backbone found"

    # Check for 17-alpha hydroxyl group pattern: C[C@H](O)
    hydroxy_17alpha_pattern = Chem.MolFromSmarts("[C;R1][C@H](O)[C;R1]")
    if not mol.HasSubstructMatch(hydroxy_17alpha_pattern):
        return False, "No 17alpha-hydroxy group found"
    
    return True, "Contains typical steroid backbone with 17alpha-hydroxy group"

# Example usage:
# print(is_17alpha_hydroxy_steroid("CC(=O)C1CCC2C1(CC(C3C2CCC4(C3CCC5C4(CCC(O)5C)(C)C=O)O)C)O"))