"""
Classifies: CHEBI:61703 nonclassic icosanoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

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
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Allow for a slightly wider range for carbon count
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (18 <= c_count <= 25):
        return False, f"Carbon count out of range: {c_count}"

    # Check for diverse oxygenation patterns
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")  # -OH
    epoxy_pattern = Chem.MolFromSmarts("[C@H]1O[C@H]1")  # Epoxy group
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")  # Carboxyl

    has_hydroxyl = mol.HasSubstructMatch(hydroxyl_pattern)
    has_epoxy = mol.HasSubstructMatch(epoxy_pattern)
    has_carboxyl = mol.HasSubstructMatch(carboxyl_pattern)

    # Ensure at least two types of oxygen-containing functional groups
    if sum([has_hydroxyl, has_epoxy, has_carboxyl]) < 2:
        return False, "Insufficient diversity in oxygenation features"

    # Exclude typical leukotriene and prostanoid structures
    leukotriene_pattern = Chem.MolFromSmarts("CCCCC(C=CC=CC=CC=CC=CC=C)C")  # Core leukotriene pattern
    prostanoid_like_pattern = Chem.MolFromSmarts("CC(C)CCC")  # Very simplified (more research is needed)

    if mol.HasSubstructMatch(leukotriene_pattern):
        return False, "Structure matches that of typical leukotrienes"

    if mol.HasSubstructMatch(prostanoid_like_pattern):
        return False, "Structure is similar to known prostanoid arrangements"
    
    # Examine long-chain fatty acid property (presence of carboxyl)
    if not has_carboxyl:
        return False, "Missing carboxyl group, typical of fatty acids"
    
    return True, "Contains characteristics of nonclassic icosanoids: C20 carbon atoms with oxygenation patterns"