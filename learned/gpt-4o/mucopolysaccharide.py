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
    uronic_acid_pattern = Chem.MolFromSmarts("[C](=O)[O][C]")  # Carboxylate ether group on a sugar residue
    if not mol.HasSubstructMatch(uronic_acid_pattern):
        return False, "No uronic acid units found"

    # Refined search: Glycosamine substructure
    glycosamine_pattern = Chem.MolFromSmarts("[C]([O,H])[C]([O,H])[N]")  # Amine on a sugar residue
    if not mol.HasSubstructMatch(glycosamine_pattern):
        return False, "No glycosamine units found"
    
    # Sulfate ester group presence
    sulfate_ester_pattern = Chem.MolFromSmarts("[O][S](=O)(=O)[O]")  # Basic sulfate ester pattern
    if not mol.HasSubstructMatch(sulfate_ester_pattern):
        return False, "Sulfate ester groups expected but not found"

    return True, "Contains uronic acids and glycosamines with esterified sulfate groups, indicating mucopolysaccharide"