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

    # Ensure 20 carbon count (specific count for C20 fatty acid backbone)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:
        return False, "Does not have at least 20 carbons, which is required for C20 fatty acids"

    # Check for at least three double bonds, enforced for icosatrienoid backbones
    dbl_bond_count = len([bond for bond in mol.GetBonds() if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE])
    if dbl_bond_count < 3:
        return False, "Insufficient double bonds (at least three required for icosatrienes)"

    # Presence and count of essential functional groups (e.g., OH, epoxy, COOH)
    hydroxyl_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[CX4][OX2H]")))
    epoxy_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[OX2R]1[CX3][CX3]1")))
    carboxylic_acid_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("C(=O)O")))
    
    if hydroxyl_count < 1 and epoxy_count < 1 and carboxylic_acid_count < 1:
        return False, "Lacking essential functional groups (need at least one OH, epoxy, or COOH group)"

    # Exclude classic leukotrienes or prostanoids
    leukotriene_exclusion = Chem.MolFromSmarts("[CH2X4]=[CHX3][CH2X4][CHX4]=[CHX3][CHX4]C") 
    prostanoid_exclusion = Chem.MolFromSmarts("C1=CCCCC1[C](=O)O")

    if mol.HasSubstructMatch(leukotriene_exclusion) or mol.HasSubstructMatch(prostanoid_exclusion):
        return False, "Matches exclusion patterns for classic icosanoids (leukotriene or prostanoid)"

    return True, "Molecule is a nonclassic icosanoid"