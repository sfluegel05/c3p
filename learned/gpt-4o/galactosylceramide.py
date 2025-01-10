"""
Classifies: CHEBI:36498 galactosylceramide
"""
from rdkit import Chem

def is_galactosylceramide(smiles: str):
    """
    Determines if a molecule is a galactosylceramide based on its SMILES string.
    A galactosylceramide is defined by a galactose head group and a long-chain
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

    # Recognize beta-D-galactose moiety pattern
    galactose_pattern = Chem.MolFromSmarts("OC[C@@H]1O[C@H](O)[C@H](O)[C@@H](CO)C1O")
    if not mol.HasSubstructMatch(galactose_pattern):
        return False, "No beta-D-galactose residue found"

    # Recognize amide linkage pattern (C(=O)N)
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide linkage found"
    
    # Detect a long-chain sphingosine backbone: [NH](C1CCCCC1)([C@H](O))
    # Typically a long carbon chain with some functionalization
    sphingosine_pattern = Chem.MolFromSmarts("*-N-*-[C@H](O)-[C@H](O)-C=C") 
    if not mol.HasSubstructMatch(sphingosine_pattern):
        return False, "Sphingosine-like backbone not detected"

    # Ensuring no phosphate groups (filter them out)
    phosphate_pattern = Chem.MolFromSmarts("OP(=O)(O)O")
    if mol.HasSubstructMatch(phosphate_pattern):
        return False, "Contains phosphate group, not a galactosylceramide"

    return True, "Contains beta-D-galactose, amide linkage, and sphingosine backbone"