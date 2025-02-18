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

    # Check for exactly two rings, as disaccharides should typically consist of two sugar rings
    if rdMolDescriptors.CalcNumRings(mol) != 2:
        return False, "Requires exactly two ring structures typical in disaccharides"
        
    # Define a glycosidic bond pattern - oxygen linking two anomeric carbons
    glycosidic_pattern = Chem.MolFromSmarts("O[C;R][C;R]")
    if not mol.HasSubstructMatch(glycosidic_pattern):
        return False, "No clear glycosidic bond pattern found"

    # Count typical atoms: carbon, oxygen; allow for substitutions
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    # Adjust carbon and oxygen limits slightly to account for common substituents in natural products
    if carbon_count < 10 or carbon_count > 30:
        return False, f"Unexpected number of carbons: {carbon_count} for disaccharide"
    if oxygen_count < 6 or oxygen_count > 18:
        return False, f"Unexpected number of oxygens: {oxygen_count} for disaccharide"

    # Check for sugar features: minimum number of hydroxyl groups for typical sugars
    if len(mol.GetSubstructMatches(Chem.MolFromSmarts("[OX2H]"))) < 6:
        return False, "Too few hydroxyl groups for a sugar-like structure"

    return True, "Contains features characteristic of disaccharides: two sugar units linked by a glycosidic bond"