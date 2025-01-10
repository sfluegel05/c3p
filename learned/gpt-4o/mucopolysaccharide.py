"""
Classifies: CHEBI:37395 mucopolysaccharide
"""
from rdkit import Chem

def is_mucopolysaccharide(smiles: str):
    """
    Determines if a molecule is a mucopolysaccharide based on its SMILES string.
    Mucopolysaccharides are polysaccharides composed of alternating units from uronic acids 
    and glycosamines, often esterified with sulfuric acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a mucopolysaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Refined search: Uronic acid substructure
    uronic_acid_pattern = Chem.MolFromSmarts("[C](=O)[O][C@,C@@H]1([O,H])[C@,C@@H]([O,H])[C@,C@@H]([O,H])[C@@H]([N]([C=O]))1")  # Pattern considering sugar ring structure
    if not mol.HasSubstructMatch(uronic_acid_pattern):
        return False, "No uronic acid units found"

    # Refined search: Glycosamine substructure
    glycosamine_pattern = Chem.MolFromSmarts("[C@,C@@H]1([O,H])[C@,C@@H]([O,H])[C@,C@@H]([O,H])[C@,C@@H]([N]([C=O]))1")  # Adjust sugar-like ring with amine
    if not mol.HasSubstructMatch(glycosamine_pattern):
        return False, "No glycosamine units found"
    
    # Sulfate ester group presence is optional; check for it
    sulfate_ester_pattern = Chem.MolFromSmarts("[O][S](=O)(=O)[O]")  # Basic sulfate ester pattern
    if not mol.HasSubstructMatch(sulfate_ester_pattern):
        message = "Contains uronic acids and glycosamines, could be a mucopolysaccharide, but lacks observed sulfate ester groups"
    else:
        message = "Contains uronic acids and glycosamines with optional sulfate ester groups, indicating mucopolysaccharide"

    return True, message