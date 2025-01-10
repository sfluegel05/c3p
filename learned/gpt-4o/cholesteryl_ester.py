"""
Classifies: CHEBI:17002 cholesteryl ester
"""
from rdkit import Chem

def is_cholesteryl_ester(smiles: str):
    """
    Determines if a molecule is a cholesteryl ester based on its SMILES string.
    A cholesteryl ester is defined by esterification of the 3-hydroxy group of cholesterol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cholesteryl ester, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Expand the cholesterol core pattern
    steroid_pattern = Chem.MolFromSmarts("C[C@H]1CCC2C(C1)CC[C@H]3[C@H]4CC[C@@H](C4)C[C@](C2)=[C@H]3")
    
    # Define expanded ester linkage pattern (flexible enough)
    ester_pattern = Chem.MolFromSmarts("C(=O)O[C@H]1CC[C@]2(C)C[C@H]3CC[C@H]4[C@H](CC[C@]4(C)C3)C2")

    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid core found representative of cholesteryl"

    ester_matches = mol.GetSubstructMatches(ester_pattern)
    steroid_matches = mol.GetSubstructMatches(steroid_pattern)
    
    # Checking if the ester linkage directly involves a steroid core hydroxyl group
    for ester_match in ester_matches:
        ester_atom_indices = set(ester_match)
        for steroid_match in steroid_matches:
            steroid_atom_indices = set(steroid_match)
            if ester_atom_indices & steroid_atom_indices: # Checks if patterns overlap at critical points
                return True, "Contains cholesterol backbone with ester functionality indicative of cholesteryl ester"
    
    return False, "Ester linkage not appropriately connected to cholesterol core"