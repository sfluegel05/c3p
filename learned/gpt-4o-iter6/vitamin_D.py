"""
Classifies: CHEBI:27300 vitamin D
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_vitamin_D(smiles: str):
    """
    Determines if a molecule is a vitamin D compound based on its SMILES string.
    Vitamin D compounds are fat-soluble hydroxy seco-steroids.
    These typically have a broken B-ring (seco-steroid) structure and specific stereochemistry.

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

    # A broader pattern capturing the typical seco-steroid structure of vitamin D.
    # This includes a triene system and common substituents.
    seco_steroid_pattern = Chem.MolFromSmarts("C1=CCC(=CC1)C=C2C[C@H](O)CC[C@]3(C3=CC=C4[C@H](CC[C@@H](O)C4)C)C2")

    if seco_steroid_pattern is None:
        return (None, None)  # Avoid failure if the pattern itself wasn't set up correctly

    if not mol.HasSubstructMatch(seco_steroid_pattern):
        return False, "Seco-steroid structure typical of vitamin D not identified"

    # Check for hydroxyl groups, which can vary but must be present.
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]") # Detect OH groups
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 1:
        return False, f"Insufficient hydroxyl groups found, at least 1 expected"

    # Verify molecular lipophilicity, modifying the threshold to avoid false negatives.
    mol_logP = rdMolDescriptors.CalcCrippenDescriptors(mol)[0]
    if mol_logP < 4:
        return False, "Molecule not sufficiently hydrophobic for vitamin D classification"

    # Check molecular weight to ensure it's within vitamin D range.
    mol_weight = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_weight < 360 or mol_weight > 600:
        return False, "Molecular weight not typical for vitamin D compounds"

    return True, "Matches expected vitamin D seco-steroid structure with hydroxyl groups and lipophilicity"