"""
Classifies: CHEBI:61703 nonclassic icosanoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nonclassic_icosanoid(smiles: str):
    """
    Determines if a molecule is a nonclassic icosanoid based on its SMILES string.
    A nonclassic icosanoid is a biologically active signaling molecule derived from 
    the oxygenation of C20 fatty acids, excluding leukotrienes and prostanoids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a nonclassic icosanoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count the number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 20:
        return False, "Does not have 20 carbons"

    # Check for presence of double bonds (C=C)
    dbl_bond_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("C=C")))
    if dbl_bond_count < 3:
        return False, "Does not have enough double bonds typical of icosanoids"

    # Look for hydroxyl (-OH), epoxy groups (O-), and carboxylic acid (-COOH) groups
    hydroxyl_pattern = Chem.MolFromSmarts("[CX4][OX2H1]")
    epoxy_pattern = Chem.MolFromSmarts("[OX2R]1[CX3][CX3]1")
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[OX2H1]")
    
    has_hydroxyl = mol.HasSubstructMatch(hydroxyl_pattern)
    has_epoxy = mol.HasSubstructMatch(epoxy_pattern)
    has_carboxylic_acid = mol.HasSubstructMatch(carboxylic_acid_pattern)

    if not (has_hydroxyl or has_epoxy or has_carboxylic_acid):
        return False, "Missing characteristic functional groups (e.g., OH, epoxy, COOH)"

    # Check for classic icosanoids exclusion
    leukotriene_pattern = Chem.MolFromSmarts("CCCCC[C@H](O)C=C(C)CCC")
    prostanoid_pattern = Chem.MolFromSmarts("CCC(=O)O")
    if mol.HasSubstructMatch(leukotriene_pattern) or mol.HasSubstructMatch(prostanoid_pattern):
        return False, "Molecule is a classic icosanoid (leukotriene or prostanoid)"

    return True, "Molecule is a nonclassic icosanoid"