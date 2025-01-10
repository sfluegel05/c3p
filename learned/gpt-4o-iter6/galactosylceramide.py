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
        return None, "Invalid SMILES string"

    # Updated sphingosine or phytosphingosine backbone patterns
    sphingosine_patterns = [
        Chem.MolFromSmarts("N[C@@H](CO)COC1O[C@H](CO)[C@H](O)[C@@H]([C@@H]1O)O"),    # Sphingosine backbone
        Chem.MolFromSmarts("N[C@H](CO)COC1O[C@@H](CO)[C@H](O)[C@H](O)[C@@H]1O"),   # Phytosphere backbone
    ]
    
    if not any(mol.HasSubstructMatch(pattern) for pattern in sphingosine_patterns):
        return False, "No sphingosine backbone found"

    # Long hydrocarbon chains attached via amide linkage
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide linkage found"

    # Flexible pattern for galactose (including variations: with or without sulfate)
    galactose_patterns = [
        Chem.MolFromSmarts("CO[C@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@@H]1O"),         # beta-D-Galactose
        Chem.MolFromSmarts("CO[C@@H]1O[C@@H](CO)[C@H](O)[C@@H](O)[C@H]1O"),       # alpha-D-Galactose
        Chem.MolFromSmarts("C1[C@H](O)[C@@H](OS(=O)(=O)O)[C@H](O)[C@@H](CO)O1")   # Sulfo beta-D-Galactose
    ]
    
    if not any(mol.HasSubstructMatch(pattern) for pattern in galactose_patterns):
        return False, "No flexible galactose head group found"

    return True, "Contains a sphingosine backbone with long-chain amide linkage and a galactose head group"