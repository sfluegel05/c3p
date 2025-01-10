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

    # Recognize galactose moiety pattern (generalized to match more cases)
    galactose_pattern = Chem.MolFromSmarts("OC[C@H]1OC(O)C(O)C(O)C1")  # Less specific to capture more variants
    if not mol.HasSubstructMatch(galactose_pattern):
        return False, "No galactose residue found"

    # Look for amide (C(=O)N) linkage pattern
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide linkage found"
    
    # Detect a long-chain backbone: generalized sphingosine detection
    sphingosine_pattern = Chem.MolFromSmarts("[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-N-[C@H](O)-C=C")  # Looking for a typical backbone
    if not mol.HasSubstructMatch(sphingosine_pattern):
        return False, "Sphingosine-like backbone not detected"

    # Ensuring absence of phosphate groups (filter them out)
    phosphate_pattern = Chem.MolFromSmarts("OP(=O)(O)O")
    if mol.HasSubstructMatch(phosphate_pattern):
        return False, "Contains phosphate group, not a galactosylceramide"
    
    # Exclude non-standard linkages and ensure sulfonate handling if necessary
    sulfonate_pattern = Chem.MolFromSmarts("S(=O)(=O)O")
    if mol.HasSubstructMatch(sulfonate_pattern):
        return False, "Contains unexpected sulfonate group"

    return True, "Contains galactose, amide linkage, and sphingosine backbone"