"""
Classifies: CHEBI:61703 nonclassic icosanoid
"""
from rdkit import Chem

def is_nonclassic_icosanoid(smiles: str):
    """
    Determines if a molecule is a nonclassic icosanoid based on its SMILES string.
    Nonclassic icosanoids are C20 fatty acid derivatives with specific oxygenation patterns,
    excluding typical leukotrienes and prostanoids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nonclassic icosanoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for approximately 20 carbon atoms (e.g., allowing carbon chain modifications)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 18 or c_count > 22:
        return False, f"Carbon count is {c_count}, typically expected around 20"

    # Check for diverse oxygenation patterns
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")  # -OH
    epoxy_pattern = Chem.MolFromSmarts("[C@H]1O[C@H]1")  # Epoxy group
    carbonyl_pattern = Chem.MolFromSmarts("[CX3](=O)")  # General carbonyl C=O group
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[OX1H0,-]")  # Carboxyl group

    has_hydroxyl = mol.HasSubstructMatch(hydroxyl_pattern)
    has_epoxy = mol.HasSubstructMatch(epoxy_pattern)
    has_carbonyl = mol.HasSubstructMatch(carbonyl_pattern)
    has_carboxyl = mol.HasSubstructMatch(carboxyl_pattern)

    # Ensure presence of at least two diverse oxygen functional groups
    if not (has_hydroxyl and has_epoxy) and not (has_carboxyl or has_carbonyl):
        return False, "Insufficient oxygenation features for a nonclassic icosanoid"

    # Exclude typical leukotriene and prostanoid structures more carefully
    supercyle_warning = Chem.MolFromSmiles('C1CCCCC1')
    leukotriene_warning = Chem.MolFromSmarts('C=CC=CC=')  # Simplified constraint for linear triene
    if mol.HasSubstructMatch(supercyle_warning) or mol.HasSubstructMatch(leukotriene_warning):
        return False, "Contains features of typical leukotriene/prostanoid structures"

    return True, "Contains characteristics of nonclassic icosanoids: C20 carbon atoms with oxygenation patterns"