"""
Classifies: CHEBI:36500 glucosylceramide
"""
from rdkit import Chem

def is_glucosylceramide(smiles: str):
    """
    Determines if a molecule is a glucosylceramide based on its SMILES string.
    A glucosylceramide is a cerebroside with a glucose head group linked to a sphingosine and a fatty acid.

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

    # Pattern for glucose head group (flexible stereochemistry)
    glucose_pattern = Chem.MolFromSmarts("O[C@H]1[C@H](O)[C@@H](O)[C@@H](CO)O[C@H]1CO")
    if not mol.HasSubstructMatch(glucose_pattern):
        return False, "No beta-D-glucosyl unit found"

    # Sphingosine backbone with amide linkage and primary alcohol
    sphingosine_pattern = Chem.MolFromSmarts("N[C@H](CO[C@H]1O[C@H](CO)C(O)[C@H](O)[C@H]1O)C(=O)")
    if not mol.HasSubstructMatch(sphingosine_pattern):
        return False, "No sphingosine structure with amide linkage found"
    
    # Detect long aliphatic chain (indicative of fatty acid)
    long_chain_pattern = Chem.MolFromSmarts("C(C(C(C(C(C(C(C(C(C(C(C"))")
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "Suitable long aliphatic chain not found"

    return True, "Contains beta-D-glucosyl unit, sphingosine backbone with amide linkage, and long fatty acid chain"