"""
Classifies: CHEBI:16158 steroid sulfate
"""
from rdkit import Chem

def is_steroid_sulfate(smiles: str):
    """
    Determines if a molecule is a steroid sulfate based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a steroid sulfate, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for sulfate groups
    sulfate_pattern = Chem.MolFromSmarts("OS(=O)(=O)O")
    
    if not mol.HasSubstructMatch(sulfate_pattern):
        return False, "No sulfate groups found"

    # Look for flexible steroid backbone pattern:
    # 3 six-membered carbon rings and 1 five-membered carbon ring. 
    # Allow optional oxygens and typical fusion points. 
    # Permit AEDQR steroid skeletons in linear arrangement with fusion at C9-C11 such as C1C[C@H]2C[C@H]3CC[C@H]4C([C@@H]4[C@@H](C3)CC2)(C1) etc.
    steroid_pattern = Chem.MolFromSmarts("C1C[C@H]2C[C@H]3CC[C@H]4[C@@H]34[C@@H](C2)CC1")
    
    if not mol.HasSubstructMatch(steroid_pattern):
        # Try another pattern in case the initial is too strict (can extend this based on observation of known errors)
        steroid_pattern_alt = Chem.MolFromSmarts("C1CCC2C3CCC4=C[C@H]([C@H]34)C2=C1") 
        if not mol.HasSubstructMatch(steroid_pattern_alt):
            return False, "No steroid backbone found (tried multiple patterns)"

    return True, "Contains steroid backbone with sulfate group(s)"