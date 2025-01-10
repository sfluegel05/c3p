"""
Classifies: CHEBI:36498 galactosylceramide
"""
from rdkit import Chem

def is_galactosylceramide(smiles: str):
    """
    Determines if a molecule is a galactosylceramide based on its SMILES string.
    A galactosylceramide is a cerebroside with a galactose monosaccharide head group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a galactosylceramide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for sphingosine-like backbone pattern
    sphingosine_pattern = Chem.MolFromSmarts("C(=C[CH:3])[CH:2](O)[CH:1](NC=O)")
    if not mol.HasSubstructMatch(sphingosine_pattern):
        return False, "No sphingosine-like backbone found"
    
    # Look for amide linkage
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide linkage found"
    
    # Look for beta-D-galactose head group (considering variability)
    galactose_pattern = Chem.MolFromSmarts("COC1OC(CO)C(O)C(O)C1O")
    if not mol.HasSubstructMatch(galactose_pattern):
        return False, "No galactose head group found"
    
    # If all the substructures are found, classify as a galactosylceramide
    return True, "Contains sphingosine-like backbone with amide linkage and galactose head group"

# Test the function with a SMILES example
smiles_example = 'C(=C/CCCCCCCCCCCCC)\\[C@@H](O)[C@@H](NC(=O)CCCCCCC/C=C\\CCCCCCCC)COC1[C@@H]([C@H]([C@H]([C@H](O1)CO)O)O)O'
result, reason = is_galactosylceramide(smiles_example)
print(f"Classification: {result}, Reason: {reason}")