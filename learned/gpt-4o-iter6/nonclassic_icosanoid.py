"""
Classifies: CHEBI:61703 nonclassic icosanoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nonclassic_icosanoid(smiles: str):
    """
    Determines if a molecule is a nonclassic icosanoid based on its SMILES string.
    These are C20 fatty acid derivatives with oxygenation, excluding leukotrienes and prostanoids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nonclassic icosanoid, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check the overall carbon count approximately equals 20
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (18 <= c_count <= 22):
        return False, f"Carbon count out of range: {c_count}"

    # Check for oxygenation: multiple types (e.g. hydroxyl, epoxy, carboxyl)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    epoxy_pattern = Chem.MolFromSmarts("[C@H]1O[C@H]1")
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    
    # At least one must match indicating oxygenation
    has_hydroxyl = mol.HasSubstructMatch(hydroxyl_pattern)
    has_epoxy = mol.HasSubstructMatch(epoxy_pattern)
    has_carboxyl = mol.HasSubstructMatch(carboxyl_pattern)

    if not (has_hydroxyl or has_epoxy or has_carboxyl):
        return False, "Insufficient oxygenation features"

    # Exclude typical leukotriene and prostanoid structures
    leukotriene_pattern = Chem.MolFromSmarts("CCCCC(C=CC=CC=CC=CC=CC=C)C")
    prostanoid_pattern = Chem.MolFromSmarts("CC(C)CCC1C2CCC3C1C(CCC3C2)C(O)=O")

    if mol.HasSubstructMatch(leukotriene_pattern):
        return False, "Structure matches that of typical leukotrienes"

    if mol.HasSubstructMatch(prostanoid_pattern):
        return False, "Structure matches that of typical prostanoids"
    
    # Check for carboxyl group indicating fatty acid
    if not has_carboxyl:
        return False, "Missing carboxyl group, typical of fatty acids"

    return True, "Contains characteristics of nonclassic icosanoids: C20 carbon atoms with oxygenation patterns"