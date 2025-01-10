"""
Classifies: CHEBI:17002 cholesteryl ester
"""
from rdkit import Chem

def is_cholesteryl_ester(smiles: str):
    """
    Determines if a molecule is a cholesteryl ester based on its SMILES string.
    A cholesteryl ester is defined by the condensation of the carboxy group of any carboxylic acid
    with the 3-hydroxy group of cholesterol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cholesteryl ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more general cholesterol core structural pattern
    # Capture fundamental sterol skeleton with relaxed stereochemical constraints
    sterol_pattern = Chem.MolFromSmarts("C1CCC2C(C1)CCC3C2CCC4(C3=CC=C4)C")
    if not mol.HasSubstructMatch(sterol_pattern):
        return False, "No cholesterol backbone found"
    
    # Ensure there is an ester linkage present
    # Looking for carbonyloxy group -OC(=O) characteristic of esters
    ester_pattern = Chem.MolFromSmarts("OC(=O)C")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester linkage found"

    # Check if the ester is attached involving the sterol core hydroxyl 
    linkage = Chem.MolFromSmarts("[C@@H]1(CC[C@@H]2C[C@@H](CC=O)CC2[C@@H]1)")
    for match in mol.GetSubstructMatches(linkage):
        # Check if any match overlaps cholesterol core, if so valid
        if any(atom_idx in match for atom_idx in mol.GetSubstructMatch(sterol_pattern)):
            return True, "Contains cholesterol backbone with ester linkage, indicating a cholesteryl ester"

    return False, "Ester linkage not appropriately connected to the cholesterol skeleton"

# Test the function with an example cholesteryl ester SMILES
smiles_example = "CC(C)CCC[C@@H](C)[C@H]1CC[C@H]2[C@@H]3CC=C4C[C@H](CC[C@]4(C)[C@H]3CC[C@]12C)OC(=O)CCCCCCC\C=C/C\C=C/CCCCC"
result = is_cholesteryl_ester(smiles_example)
print(result)