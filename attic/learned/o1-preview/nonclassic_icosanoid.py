"""
Classifies: CHEBI:61703 nonclassic icosanoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nonclassic_icosanoid(smiles: str):
    """
    Determines if a molecule is a nonclassic icosanoid based on its SMILES string.
    A nonclassic icosanoid is any biologically active signalling molecule made by oxygenation of C20 fatty acids 
    other than the classic icosanoids (the leukotrienes and the prostanoids).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nonclassic icosanoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count the number of carbon atoms
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons < 18 or num_carbons > 36:
        return False, f"Molecule has {num_carbons} carbons, expected between 18 and 36"

    # Check for carboxylic acid group (-COOH) or ester group (-COOR)
    carboxylic_acid_or_ester_pattern = Chem.MolFromSmarts('C(=O)[O;H1,H0]')
    if not mol.HasSubstructMatch(carboxylic_acid_or_ester_pattern):
        return False, "No carboxylic acid or ester group found"

    # Check for oxygenated functional groups (hydroxyl, epoxy, ketone, peroxide, ether)
    oxygen_functional_groups = [
        Chem.MolFromSmarts('[OX2H]'),               # hydroxyl group
        Chem.MolFromSmarts('[CX3](=O)[OX2H0][CX3](=O)'),  # peroxide group
        Chem.MolFromSmarts('C1OC1'),                # epoxide ring
        Chem.MolFromSmarts('C=O'),                  # ketone or aldehyde
        Chem.MolFromSmarts('[CX4][OX2][CX4]'),      # ether linkage
    ]

    oxygenated = False
    for pattern in oxygen_functional_groups:
        if mol.HasSubstructMatch(pattern):
            oxygenated = True
            break
    if not oxygenated:
        return False, "No oxygenated functional groups found"

    # Check for multiple double bonds (more than 2)
    num_double_bonds = len([bond for bond in mol.GetBonds() if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE])
    if num_double_bonds < 3:
        return False, f"Only {num_double_bonds} double bonds found, expected multiple double bonds"

    # Exclude classic icosanoids (prostanoids and leukotrienes)
    # Prostanoids have a cyclopentane ring with adjacent oxygenated groups
    prostanoid_pattern = Chem.MolFromSmarts('C1CCC(C1)C(=O)')  # Simplified prostanoid ring with ketone
    if mol.HasSubstructMatch(prostanoid_pattern):
        return False, "Molecule matches prostanoid pattern (classic icosanoid)"

    # Leukotrienes have a conjugated triene system and terminal carboxylic acid
    leukotriene_pattern = Chem.MolFromSmarts('C=CC=CC=CC(=O)O')  # Simplified leukotriene pattern
    if mol.HasSubstructMatch(leukotriene_pattern):
        return False, "Molecule matches leukotriene pattern (classic icosanoid)"

    return True, "Molecule is a nonclassic icosanoid"