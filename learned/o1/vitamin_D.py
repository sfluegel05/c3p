"""
Classifies: CHEBI:27300 vitamin D
"""
"""
Classifies: Vitamin D compounds
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_vitamin_D(smiles: str):
    """
    Determines if a molecule is a vitamin D compound based on its SMILES string.
    Vitamin D compounds are secosteroids with a broken B-ring, a conjugated triene system,
    hydroxyl groups at specific positions, and a hydrocarbon side chain at position 17.

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

    # Check for steroid backbone (four fused rings)
    steroid_smarts = 'C1CCC2C(C1)CCC3C2CCC4=C3C=CC=C4'  # Simplified steroid nucleus
    steroid_pattern = Chem.MolFromSmarts(steroid_smarts)
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"

    # Check for secosteroid feature (broken B-ring between carbons 9 and 10)
    # Looking for a steroid backbone missing one ring (secosteroid)
    secosteroid_smarts = 'C1CCC2C(C1)CC=C3C2CCC4=C3C=CC=C4'  # Opened B-ring
    secosteroid_pattern = Chem.MolFromSmarts(secosteroid_smarts)
    if not mol.HasSubstructMatch(secosteroid_pattern):
        return False, "No secosteroid core (opened B-ring) found"

    # Check for conjugated triene system
    triene_smarts = 'C=C-C=C-C=C'  # Conjugated triene
    triene_pattern = Chem.MolFromSmarts(triene_smarts)
    if not mol.HasSubstructMatch(triene_pattern):
        return False, "No conjugated triene system found"

    # Check for at least two hydroxyl groups
    hydroxyl_smarts = '[OX2H]'  # Hydroxyl group
    hydroxyl_pattern = Chem.MolFromSmarts(hydroxyl_smarts)
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 2:
        return False, f"Found {len(hydroxyl_matches)} hydroxyl group(s), need at least 2"

    # Estimate molecular weight for a vitamin D range (approximately 380-420 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 380 or mol_wt > 450:
        return False, f"Molecular weight {mol_wt:.2f} Da not in typical range for vitamin D"

    # Check for side chain at position 17 (D-ring)
    # Looking for a hydrocarbon chain attached to the steroid nucleus
    side_chain_smarts = 'C[C@H](CCCC(C)C)C3CCC4C2C(CCC4)=C1C=CC=C1C2C3'  # Side chain pattern
    side_chain_pattern = Chem.MolFromSmarts(side_chain_smarts)
    if not mol.HasSubstructMatch(side_chain_pattern):
        return False, "No appropriate side chain found at position 17"

    return True, "Molecule matches key structural features of vitamin D"