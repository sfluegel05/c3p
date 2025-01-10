"""
Classifies: CHEBI:27300 vitamin D
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_vitamin_D(smiles: str):
    """
    Determines if a molecule is a vitamin D compound based on its SMILES string.
    Vitamin D compounds are defined as fat-soluble hydroxy seco-steroids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a vitamin D compound, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define an improved seco-steroid pattern specific to vitamin D
    # Key features are the rings A,B, and a break in C/D rings for vitamin D
    seco_steroid_pattern = Chem.MolFromSmarts("C1CCC2[C@@H]3CC[C@H](O)C(C)=C3C=C2C1")  # Representative pattern

    if seco_steroid_pattern is None:
        return (None, None)  # Avoid failure if the pattern itself wasn't set up correctly

    if not mol.HasSubstructMatch(seco_steroid_pattern):
        return False, "Seco-steroid backbone typical of vitamin D not identified"

    # Check for presence of hydroxyl groups (-OH)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")  # OH group
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 1:
        return False, f"Insufficient hydroxyl groups found, at least 1 expected"

    # Adjusted check to ensure the molecule is hydrophobic/lipophilic enough
    mol_logP = rdMolDescriptors.CalcCrippenDescriptors(mol)[0]
    if mol_logP < 4:  # Adjusted threshold
        return False, "Molecule not sufficiently hydrophobic for vitamin D classification"

    # Check molecular weight to further filter likely vitamin D compounds
    mol_weight = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_weight < 350:
        return False, "Molecular weight lower than typical vitamin D compounds"

    return True, "Matches expected vitamin D seco-steroid structure with hydroxyl groups and lipophilicity"