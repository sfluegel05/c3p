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

    # Primary pattern for identifying the core seco-steroid structure
    # Vitamin D compounds typically have the A- and C-rings intact with a Seco-B-ring structure, sometimes a triene.
    seco_B_ring_pattern = Chem.MolFromSmarts("C1CCC2C(C1)CCC=C2C=C")
    if seco_B_ring_pattern is None:
        return (None, None)  # Avoid failure if the pattern itself wasn't set up correctly

    # Match for the primary structure
    matched_primary = mol.HasSubstructMatch(seco_B_ring_pattern)

    # Check for hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 2:  # At least two hydroxy groups are necessary for primary vitamin D
        return False, "Insufficient hydroxyl groups found, at least 2 expected"

    # Confirm hydrophobic characteristics
    mol_logP = rdMolDescriptors.CalcCrippenDescriptors(mol)[0]
    if mol_logP < 4:
        return False, "Molecule not sufficiently hydrophobic for vitamin D classification"

    # Confirm molecular weight
    mol_weight = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_weight < 360 or mol_weight > 580:
        return False, "Molecular weight not typical for vitamin D compounds"

    if matched_primary:
        return True, "Matches typical vitamin D seco-steroid structure"

    # If primary matching fails, return non-classification
    return False, "Does not match vitamin D chemical description"