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

    # Check for at least two rings, allowing for some flexibility in complex disaccharides
    if rdMolDescriptors.CalcNumRings(mol) < 2:
        return False, "Too few ring structures to represent a disaccharide"
        
    # Define a pattern indicating glycosidic linkage (oxygen linking two sugar-like units)
    glycosidic_pattern = Chem.MolFromSmarts("O[C;R][C;R]")
    if not mol.HasSubstructMatch(glycosidic_pattern):
        return False, "No clear glycosidic bond pattern found"

    # Count typical atoms for disaccharides: carbon, oxygen, optional nitrogen (for derivatives)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    nitrogen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)

    # Loosen the constraints, given some derivatives or substitutions: note commonly with 10-24 carbons and 6-12 oxygens
    if carbon_count < 10 or carbon_count > 28:
        return False, f"Unexpected number of carbons: {carbon_count} for disaccharide"
    if oxygen_count < 6 or oxygen_count > 16:
        return False, f"Unexpected number of oxygens: {oxygen_count} for disaccharide"

    # Check for sugars-like features, such as hydroxyl groups or ethers
    if len(mol.GetSubstructMatches(Chem.MolFromSmarts("[OX2H]"))) < 4:
        return False, "Too few hydroxyl groups for a sugar-like structure"
    if not mol.HasSubstructMatch(Chem.MolFromSmarts("[OX2][CX4]")):
        return False, "Lack of typical ether linkage in a disaccharide"

    return True, "Contains features characteristic of disaccharides: two sugar units linked by a glycosidic bond"