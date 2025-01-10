"""
Classifies: CHEBI:36498 galactosylceramide
"""
from rdkit import Chem

def is_galactosylceramide(smiles: str):
    """
    Determines if a molecule is a galactosylceramide based on its SMILES string.
    A galactosylceramide is a cerebroside with a galactose monosaccharide head group.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is a galactosylceramide, False otherwise.
        str: Reason for classification.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # General sphingosine backbone pattern
    sphingosine_pattern = Chem.MolFromSmarts("C[C@H](O)CC(N)[C@H](O)")
    if not mol.HasSubstructMatch(sphingosine_pattern):
        return False, "No general sphingosine backbone found"

    # Amide linkage - long hydrocarbon chain
    amide_pattern = Chem.MolFromSmarts("C(=O)N[C@@H](C)CCCCCCCCCCCCCCCC")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No long-chain amide linkage found"

    # Flexible pattern for galactose including variations like sulfation
    galactose_pattern = Chem.MolFromSmarts("CO[C@H]1O[C@H](CO)C(O)C(O)C1O")
    if not mol.HasSubstructMatch(galactose_pattern):
        return False, "No flexible galactose head group found"

    return True, "Contains a general sphingosine backbone with long-chain amide linkage and a flexible galactose head group"