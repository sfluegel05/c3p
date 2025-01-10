"""
Classifies: CHEBI:35915 sterol ester
"""
from rdkit import Chem
from rdkit.Chem import rdChemReactions

def is_sterol_ester(smiles: str):
    """
    Determines if a molecule is a sterol ester based on its SMILES string.
    A sterol ester is a steroid ester obtained by formal condensation of 
    the carboxy group of any carboxylic acid with the 3-hydroxy group of a sterol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sterol ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more flexible sterol backbone pattern
    # Sterol: Cyclopentanophenanthrene with oxygens on the ring
    sterol_pattern = Chem.MolFromSmarts("C1CCC2(C(C1)CC3C(C2)CC4C(CC(CC4)C3)O)C")
    if not mol.HasSubstructMatch(sterol_pattern):
        return False, "No sterol backbone found"

    # Ester linkage pattern specifically with linkage to a sterol -OH group
    # General ester: O=C-O, but linked to a C-OH in the sterol
    ester_linker_pattern = Chem.MolFromSmarts("[C;R][CX4H](O[C;R])C(=O)O")
    if not mol.HasSubstructMatch(ester_linker_pattern):
        return False, "No appropriate ester linkage found"

    return True, "Contains sterol backbone with appropriate ester linkage"

# Test the function with an example
example_smiles = "C1[C@@]2([C@]3(CC[C@]4([C@]([C@@]3(CC=C2C[C@H](C1)OC(=O)CCCCCCCCC/C=C\\C/C=C\\C/C=C\\CC)[H])(CC[C@@]4([C@H](C)CCCC(C)C)[H])[H])C)[H])C"
result, reason = is_sterol_ester(example_smiles)
print(f"Result: {result}, Reason: {reason}")