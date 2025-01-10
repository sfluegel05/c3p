"""
Classifies: CHEBI:36498 galactosylceramide
"""
from rdkit import Chem

def is_galactosylceramide(smiles: str):
    """
    Determines if a molecule is a galactosylceramide based on its SMILES string.
    A galactosylceramide includes a galactose head group and a long-chain
    ceramide backbone.

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

    # Recognize galactose moiety pattern, including sulfated versions
    galactose_pattern = Chem.MolFromSmarts("OC[C@@H]1O[C@H](CO)[C@H](O)[C@@H](O)[C@H]1O")
    if not mol.HasSubstructMatch(galactose_pattern):
        return False, "No galactose residue found"

    # Look for amide (C(=O)N) linkage pattern
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide linkage found"
    
    # Detect a sphingosine-like backbone: flexible long-chain pattern
    sphingosine_pattern = Chem.MolFromSmarts("[#6]~[CH2]~[CH2]~[CH2]~[CH]~[N]~[CH](O)~[#6]")  
    if not mol.HasSubstructMatch(sphingosine_pattern):
        return False, "Sphingosine-like backbone not detected"
    
    # Exclude components typically non-galactosylceramide except sulfonate
    phosphate_pattern = Chem.MolFromSmarts("OP(=O)(O)O")
    if mol.HasSubstructMatch(phosphate_pattern):
        return False, "Contains phosphate group, not a galactosylceramide"

    return True, "Contains galactose, amide linkage, and sphingosine backbone"