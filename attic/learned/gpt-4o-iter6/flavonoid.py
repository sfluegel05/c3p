"""
Classifies: CHEBI:47916 flavonoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_flavonoid(smiles: str):
    """
    Determines if a molecule is a flavonoid based on its SMILES string.
    A flavonoid is characterized by a 1-benzopyran structure with an aryl group at position 2.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for 1-benzopyran ring system (simplified chromene core) with an aryl at position 2
    chromene_pattern = Chem.MolFromSmarts("C1=CC=C2C(=C1)OC=C2")  # Basic 1-benzopyran
    if not mol.HasSubstructMatch(chromene_pattern):
        return False, "No 1-benzopyran core found"

    # Look for aryl substitution at position 2
    aryl_pattern = Chem.MolFromSmarts("c1ccc(cc1)-C2=CC=CO2")  # Aryl at position 2 of chromene
    if not mol.HasSubstructMatch(aryl_pattern):
        return False, "No aryl group at correct position"

    # Count hydroxyl groups
    hydroxyl_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() == 1)
    if hydroxyl_count < 1:
        return False, "Insufficient hydroxyl groups for a typical flavonoid structure"

    return True, "Contains flavonoid core with aryl substitution at position 2"