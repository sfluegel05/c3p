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
    if num_carbons < 18 or num_carbons > 22:
        return False, f"Molecule has {num_carbons} carbons, expected around 20"

    # Check for carboxylic acid group (-COOH)
    carboxylic_acid_pattern = Chem.MolFromSmarts('C(=O)[O;H1]')
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Check for oxygenated functional groups (hydroxyl, epoxy, ketone)
    oxygen_functional_groups = [
        Chem.MolFromSmarts('[OX2H]'),          # hydroxyl group
        Chem.MolFromSmarts('C1OC1'),           # epoxide ring
        Chem.MolFromSmarts('C=O'),             # ketone or aldehyde
        Chem.MolFromSmarts('[#6]-[O]-[#6]'),   # ether linkage
    ]

    oxygenated = False
    for pattern in oxygen_functional_groups:
        if mol.HasSubstructMatch(pattern):
            oxygenated = True
            break
    if not oxygenated:
        return False, "No oxygenated functional groups found"

    # Check for multiple double bonds
    num_double_bonds = len([bond for bond in mol.GetBonds() if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE])
    if num_double_bonds < 2:
        return False, f"Only {num_double_bonds} double bonds found, expected multiple double bonds"

    # Exclude classic icosanoids (prostanoids and leukotrienes)
    # Simplified patterns for prostanoids and leukotrienes
    prostanoid_pattern = Chem.MolFromSmarts('C1(=O)CC[C@H](O)CC1')  # Simplified prostanoid ring
    leukotriene_pattern = Chem.MolFromSmarts('CCCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC(=O)O')  # Simplified leukotriene chain

    if mol.HasSubstructMatch(prostanoid_pattern):
        return False, "Molecule matches prostanoid pattern (classic icosanoid)"
    if mol.HasSubstructMatch(leukotriene_pattern):
        return False, "Molecule matches leukotriene pattern (classic icosanoid)"

    return True, "Molecule is a nonclassic icosanoid"