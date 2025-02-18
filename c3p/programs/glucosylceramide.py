"""
Classifies: CHEBI:36500 glucosylceramide
"""
from rdkit import Chem

def is_glucosylceramide(smiles: str):
    """
    Determines if a molecule is a glucosylceramide based on its SMILES string.
    A glucosylceramide consists of a glucose moiety linked to a sphingosine backbone
    with a fatty acid chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glucosylceramide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify -glucose ring with stereochemistry
    glucose_pattern = Chem.MolFromSmarts("O[C@@H]1[C@H](O)[C@@H](O)[C@H](CO)O[C@H]1CO")
    if not mol.HasSubstructMatch(glucose_pattern):
        return False, "No -glucose group found"

    # Identify amide linkage (N-C(=O))
    amide_pattern = Chem.MolFromSmarts("NC(=O)")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide linkage found"

    # Identify sphingosine backbone (assumption: long hydrocarbon with amine and hydroxyl)
    sphingosine_pattern = Chem.MolFromSmarts("[C@@H](NC=O)[C@@H](O)CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O")
    if not mol.HasSubstructMatch(sphingosine_pattern):
        return False, "No sphingosine backbone found"

    # Identify long fatty acid chain (10+ carbons)
    long_chain_pattern = Chem.MolFromSmarts("CCCCCCCCCCCC")
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "No long fatty acid chain found"

    return True, "Contains glucose moiety linked to sphingosine with a fatty acid chain"