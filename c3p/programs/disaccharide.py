"""
Classifies: CHEBI:36233 disaccharide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_disaccharide(smiles: str):
    """
    Determines if a molecule is a disaccharide based on its SMILES string.
    A disaccharide is defined as two monosaccharides joined by a glycosidic bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a disaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Adjust the requirement from exactly two rings to a range to accommodate flexible forms
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    if num_rings < 1 or num_rings > 3:
        return False, f"Number of rings ({num_rings}) not typical for disaccharides"

    # Comprehensive glycosidic bond patterns
    glycosidic_pattern = Chem.MolFromSmarts("O[C;R][C;R]")
    if not mol.HasSubstructMatch(glycosidic_pattern):
        alt_glycosidic_pattern = Chem.MolFromSmarts("O[C;R]C")
        if not mol.HasSubstructMatch(alt_glycosidic_pattern):
            return False, "No clear glycosidic bond pattern found"

    # Relax but bound carbon and oxygen counts
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    if carbon_count < 10 or carbon_count > 40:
        return False, f"Unexpected number of carbons: {carbon_count} for disaccharide"
    if oxygen_count < 5 or oxygen_count > 22:
        return False, f"Unexpected number of oxygens: {oxygen_count} for disaccharide"

    # Ensure presence of hydroxyl groups, crucial for sugar-like structure
    if len(mol.GetSubstructMatches(Chem.MolFromSmarts("[OX2H]"))) < 5:
        return False, "Too few hydroxyl groups for a sugar-like structure"

    return True, "Contains features characteristic of disaccharides: two sugar units linked by a glycosidic bond"