"""
Classifies: CHEBI:33303 chalcogen
"""
"""
Classifies: Chalcogen (Group 16 elements)
"""
from rdkit import Chem

def is_chalcogen(smiles: str):
    """
    Determines if a SMILES represents a chalcogen atom (Group 16 elements: O, S, Se, Te, Po, Lv).

    Args:
        smiles (str): SMILES string

    Returns:
        bool: True if it's a chalcogen atom, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Check if the molecule consists of a single atom
    if mol.GetNumAtoms() != 1:
        return False, "Not a single atom"
    
    atom = mol.GetAtomWithIdx(0)
    atomic_num = atom.GetAtomicNum()
    
    # Group 16 elements: O(8), S(16), Se(34), Te(52), Po(84), Lv(116)
    if atomic_num in {8, 16, 34, 52, 84, 116}:
        return True, f"Element {atom.GetSymbol()} is a Group 16 chalcogen"
    else:
        return False, f"Element {atom.GetSymbol()} is not a Group 16 chalcogen"