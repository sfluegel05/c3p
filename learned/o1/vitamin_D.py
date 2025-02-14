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
    Vitamin D compounds are secosteroids with an open B-ring, fused rings C and D,
    a conjugated triene system in the A-ring, hydroxyl groups at specific positions, 
    and a hydrocarbon side chain.

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

    # Define SMARTS patterns for vitamin D features

    # Secosteroid core with open B-ring and fused rings C and D
    secosteroid_core_smarts = """
    [
        C;R2;!R3]1
        [C;R2;!R3][C;R2;!R3][C;R2;!R3][C;R2;!R3][C;R2;!R3]1
        [C;R2;!R3]
    ]"""
    secosteroid_core_pattern = Chem.MolFromSmarts("""
    [#6]1-[#6]-[#6]-[#6]-[#6]-1-[#6]
    """)

    # Open B-ring (missing in vitamin D)
    B_ring_smarts = 'C1CCC1'  # Cyclobutane ring (should not be present)
    B_ring_pattern = Chem.MolFromSmarts(B_ring_smarts)
    if mol.HasSubstructMatch(B_ring_pattern):
        return False, "B-ring (cyclobutane) is present, not a secosteroid"

    # Fused rings C (cyclehexane) and D (cyclopentane)
    fused_ring_pattern = Chem.MolFromSmarts('[#6]1-[#6]-[#6]-[#6]-[#6]-[#6]-1-[#6]2-[#6]-[#6]-[#6]-[#6]-2')
    if not mol.HasSubstructMatch(fused_ring_pattern):
        return False, "Fused rings C and D not found"

    # Conjugated triene system in the A-ring
    triene_smarts = 'C=C-C=C-C=C'  # Conjugated triene system
    triene_pattern = Chem.MolFromSmarts(triene_smarts)
    if not mol.HasSubstructMatch(triene_pattern):
        return False, "No conjugated triene system found in the A-ring"

    # Check for hydroxyl groups
    hydroxyl_smarts = '[OX2H]'  # Hydroxyl group
    hydroxyl_pattern = Chem.MolFromSmarts(hydroxyl_smarts)
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    hydroxyl_count = len(hydroxyl_matches)
    if hydroxyl_count < 2:
        return False, f"Found {hydroxyl_count} hydroxyl group(s), need at least 2"

    # Check for hydrocarbon side chain (aliphatic chain attached to ring D)
    side_chain_smarts = '[C]-[CH2]-[C](C)(C)O'  # Simplified side chain pattern
    side_chain_pattern = Chem.MolFromSmarts(side_chain_smarts)
    if not mol.HasSubstructMatch(side_chain_pattern):
        return False, "Hydrocarbon side chain not found"

    # Estimate molecular weight for vitamin D range (approximately 380-450 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 350 or mol_wt > 500:
        return False, f"Molecular weight {mol_wt:.2f} Da not in typical range for vitamin D"

    return True, "Molecule matches key structural features of vitamin D"