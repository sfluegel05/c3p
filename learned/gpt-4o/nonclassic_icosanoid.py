"""
Classifies: CHEBI:61703 nonclassic icosanoid
"""
from rdkit import Chem

def is_nonclassic_icosanoid(smiles: str):
    """
    Determines if a molecule is a nonclassic icosanoid based on its SMILES string.
    A nonclassic icosanoid is a biologically active signaling molecule derived from 
    the oxygenation of C20 fatty acids, excluding classic leukotrienes and prostanoids.

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

    # Ensure 20 carbon count (if length of carbon chain specifies the base icosanoid)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 20:
        return False, "Does not have 20 carbons"

    # Double bonds presence is crucial for the icosa structure
    dbl_bond_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("C=C")))
    if dbl_bond_count < 3:
        return False, "Insufficient double bonds"

    # Presence of functional groups (OH, epoxy, and COOH) typical in nonclassic icosanoids
    hydroxyl_pattern = Chem.MolFromSmarts("[CX4][OX2H1]")
    epoxy_pattern = Chem.MolFromSmarts("[OX2R]1[CX3][CX3]1")
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[OX2H1]")

    has_hydroxyl = mol.HasSubstructMatch(hydroxyl_pattern)
    has_epoxy = mol.HasSubstructMatch(epoxy_pattern)
    has_carboxylic_acid = mol.HasSubstructMatch(carboxylic_acid_pattern)
    
    if not (has_hydroxyl or has_epoxy or has_carboxylic_acid):
        return False, "Lacking essential functional groups (e.g., OH, epoxy, COOH)"

    # Refined exclusion for classic icosanoids (leukotrienes and prostanoids)
    # Adapted by synthesizing heuristic-based observation of known classic structures
    leukotriene_exclusion = Chem.MolFromSmarts("[CH2X4][CHX4]=[CH][CH2X4]CCCC(C)C") 
    prostanoid_exclusion = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")

    if mol.HasSubstructMatch(leukotriene_exclusion) or mol.HasSubstructMatch(prostanoid_exclusion):
        return False, "Structure aligns with classic icosanoid templates (leukotriene or prostanoid)"

    return True, "Molecule is a nonclassic icosanoid"