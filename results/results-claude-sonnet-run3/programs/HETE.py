from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.AllChem import GetBestRMS

def is_HETE(smiles: str):
    """
    Determines if a molecule is a HETE (hydroxyeicosatetraenoic acid).
    HETEs are monohydroxylated icosanoids with 4 double bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a HETE, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Check carbon count (should be 20 for eicosanoid)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')
    if carbon_count != 20:
        return False, f"Not an eicosanoid - has {carbon_count} carbons instead of 20"

    # Check for carboxylic acid group
    carboxylic_pattern = Chem.MolFromSmarts('C(=O)O')
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "Missing carboxylic acid group"

    # Count hydroxyl groups (should be exactly one, not counting the carboxylic acid)
    hydroxyl_pattern = Chem.MolFromSmarts('[CH2,CH]-[OH]')
    hydroxyl_matches = len(mol.GetSubstructMatches(hydroxyl_pattern))
    if hydroxyl_matches != 1:
        return False, f"Should have exactly one hydroxyl group, found {hydroxyl_matches}"

    # Count double bonds (should be 4)
    double_bond_pattern = Chem.MolFromSmarts('C=C')
    double_bond_matches = len(mol.GetSubstructMatches(double_bond_pattern))
    if double_bond_matches != 4:
        return False, f"Should have exactly 4 double bonds, found {double_bond_matches}"

    return True, "Molecule is a HETE (hydroxyeicosatetraenoic acid)"
# Pr=1.0
# Recall=1.0