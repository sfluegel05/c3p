"""
Classifies: CHEBI:36498 galactosylceramide
"""
from rdkit import Chem

def is_galactosylceramide(smiles: str):
    """
    Determines if a molecule is a galactosylceramide based on its SMILES string.
    Galactosylceramides are cerebrosides with a galactose head group and a sphingolipid backbone.

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

    # Detect the galactose portion, including common modifications like sulfate
    galactose_pattern = Chem.MolFromSmarts("OC[C@H]1O[C@H](O)[C@@H](O)[C@@H](~[CH2S](=O)(=O)[O])[C@H]1O | C[C@H]1O[C@H](O)[C@@H](O)[C@@H](O)[C@H]1O")
    if not mol.HasSubstructMatch(galactose_pattern):
        return False, "No recognizable galactose head group found"

    # Ensure an amide bond connecting a fatty acid to a sphingolipid
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide linkage found"
    
    # Broader sphingolipid backbone pattern with allowance for common variations
    sphingolipid_pattern = Chem.MolFromSmarts("C(=O)N[C@@H](CO[*])[CX3]=C[*]")  # Broader pattern considering sphingoid variations
    if not mol.HasSubstructMatch(sphingolipid_pattern):
        return False, "No sphingolipid backbone found"
    
    # Checking elemental composition for carbohydrates and lipid chain validity
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    if o_count < 6:
        return False, "Insufficient oxygen atoms for necessary sugar and amide functionalities"
    if c_count < 30:
        return False, "Insufficient carbon atoms for long lipid and sugar chains"

    return True, "Contains galactose head group, sphingolipid backbone, and necessary structural features"