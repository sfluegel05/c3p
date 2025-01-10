"""
Classifies: CHEBI:36498 galactosylceramide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

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
    
    # Look for sphingosine backbone pattern
    # This is typically represented as a long chain amino alcohol
    sphingosine_pattern = Chem.MolFromSmarts("C[C@@H](O)[C@H](N)CCO")
    if not mol.HasSubstructMatch(sphingosine_pattern):
        return False, "No sphingosine backbone found"
    
    # Look for amide linkage
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide linkage found"
    
    # Look for galactose head group
    # Galactose is a specific sugar which can be recognized by its stereochemistry - a hexose sugar with specific hydroxyl group orientations
    galactose_pattern = Chem.MolFromSmarts("OC[C@@H]1O[C@@H]([C@H](O)[C@H](O)[C@@H]1O)CO")
    if not mol.HasSubstructMatch(galactose_pattern):
        return False, "No galactose head group found"
    
    # If all the substructures are found, classify as a galactosylceramide
    return True, "Contains sphingosine backbone with amide linkage and galactose head group"

# Example Usage
smiles_example = 'C(=C/CCCCCCCCCCCCC)\\[C@@H](O)[C@@H](NC(=O)CCCCCCC/C=C\\CCCCCCCC)COC1[C@@H]([C@H]([C@H]([C@H](O1)CO)O)O)O'
result, reason = is_galactosylceramide(smiles_example)
print(f"Classification: {result}, Reason: {reason}")