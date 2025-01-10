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

    # Define a cholesterol backbone pattern with potential understandable stereochemistry
    # Capturing basic sterol structure with essential rings and alkyl side chain
    # Included stereochemistry flexibility by removing specific chiral specs from the smart pattern.
    cholesterol_pattern = Chem.MolFromSmarts("C1CCC2C(C1)CCC3C2CCC4(C3=CC=C4)C")
    sterol_matches = mol.GetSubstructMatches(cholesterol_pattern)
    if not sterol_matches:
        return False, "No cholesterol backbone found"
    
    # Check there is an ester group attached (-OC(=O)-)
    ester_pattern = Chem.MolFromSmarts("OC(=O)C")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester linkage found"

    # Check connection: ester linkage attached to cholesterol core (position 3 OH)
    # Find correct connection to ensure ester is linked to the steroid nucleus
    for match in sterol_matches:
        # Check if ester connects to cholesterol hydroxyl position
        connection_pattern = Chem.MolFromSmarts("C1(O[C]=O)CCC2C(C1)CCC3C2CCC4(C3=CC=C4)C")
        if mol.HasSubstructMatch(connection_pattern):
            return True, "Contains cholesterol backbone with ester linkage, indicating a cholesteryl ester"
    
    # If nothing matches return false
    return False, "Ester linkage not appropriately connected to the cholesterol skeleton"

# Example test
smiles_example = "CC(C)CCC[C@@H](C)[C@H]1CC[C@H]2[C@@H]3CC=C4C[C@H](CC[C@]4(C)[C@H]3CC[C@]12C)OC(=O)CCCCCCC\C=C/C\C=C/CCCCC"
result = is_cholesteryl_ester(smiles_example)
print(result)