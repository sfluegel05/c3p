"""
Classifies: CHEBI:18154 polysaccharide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polysaccharide(smiles):
    """
    Determines if a molecule is a polysaccharide based on its SMILES string.
    A polysaccharide consists of more than ten monosaccharide residues linked glycosidically.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polysaccharide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES string"

    # Define patterns for monosaccharide units and glycosidic linkages
    monosaccharide_pattern = Chem.MolFromSmarts("[C@@H]([C@H]([C@H](O)C=O)O)[C@@H]")  # Basic monosaccharide structure
    glycosidic_link_pattern = Chem.MolFromSmarts("O[C@H]")  # Simplified glycosidic bond
    
    # Count monosaccharide units
    mono_count = len(mol.GetSubstructMatches(monosaccharide_pattern))
    if mono_count < 10:
        return False, f"Contains {mono_count} monosaccharides, less than required 10 for polysaccharide"

    # Count glycosidic bonds
    glyco_count = len(mol.GetSubstructMatches(glycosidic_link_pattern))
    if glyco_count < mono_count - 1:
        return False, f"Insufficient glycosidic linkages ({glyco_count}) for {mono_count} monosaccharides"

    # Final classification
    return True, "Contains sufficient monosaccharide moieties linked glycosidically to classify as a polysaccharide"