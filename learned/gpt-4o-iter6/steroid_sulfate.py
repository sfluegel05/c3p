"""
Classifies: CHEBI:16158 steroid sulfate
"""
from rdkit import Chem

def is_steroid_sulfate(smiles: str):
    """
    Determines if a molecule is a steroid sulfate based on its SMILES string.
    A steroid sulfate is a sulfuric ester obtained by the condensation of a hydroxy group of any steroid with sulfuric acid.

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

    # Improved SMARTS pattern for general steroid core - more flexible structure
    steroid_pattern = Chem.MolFromSmarts("C1CCC2C(C1)CCC3C2CCC4(C3C=CC=C4)C") # Simplest common backbone

    # Check if there's a pattern matching a steroid backbone
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
      
    # SMARTS pattern to identify sulfate group bound through an ester linkage
    sulfate_pattern = Chem.MolFromSmarts("OS(=O)(=O)[O,C]-")

    # Check for sulfate linkage to hydroxy group
    if not mol.HasSubstructMatch(sulfate_pattern):
        return False, "No sulfate group properly linked"

    # If both the steroid backbone and sulfate ester linkage are found
    return True, "Contains steroid backbone with sulfate ester linkage"

# Example for testing
smiles = "[Na+].[H][C@]12CC[C@]3(C)C(=O)CC[C@@]3([H])[C@]1([H])CCc1cc(OS([O-])(=O)=O)ccc21"
result, reason = is_steroid_sulfate(smiles)
print(f"Result: {result}, Reason: {reason}")