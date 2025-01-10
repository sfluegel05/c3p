"""
Classifies: CHEBI:27325 xanthophyll
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_xanthophyll(smiles: str):
    """
    Determines if a molecule is a xanthophyll based on its SMILES string.
    Xanthophylls are oxygenated carotenes, containing oxygen atoms in addition to the typical carotenoid backbone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a xanthophyll, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of oxygen atoms
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count == 0:
        return False, "No oxygen atoms found, not a xanthophyll"

    # Check for the presence of a long polyene chain (carotenoid backbone)
    # Use a more flexible pattern to identify the carotenoid backbone
    carotenoid_pattern = Chem.MolFromSmarts("[CH2,CH,CH0]~[CH,CH0]~[CH,CH0]~[CH,CH0]~[CH,CH0]~[CH,CH0]~[CH,CH0]~[CH,CH0]")
    if not mol.HasSubstructMatch(carotenoid_pattern):
        return False, "No carotenoid backbone found"

    # Check molecular weight (xanthophylls are typically >400 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, "Molecular weight too low for a xanthophyll"

    # Check for specific oxygen-containing functional groups (hydroxyl, epoxy, carbonyl, ether, carboxyl)
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    epoxy_pattern = Chem.MolFromSmarts("[OX2]C1CC1")
    carbonyl_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")
    ether_pattern = Chem.MolFromSmarts("[OX2]C")
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[OH]")

    has_hydroxyl = mol.HasSubstructMatch(hydroxyl_pattern)
    has_epoxy = mol.HasSubstructMatch(epoxy_pattern)
    has_carbonyl = mol.HasSubstructMatch(carbonyl_pattern)
    has_ether = mol.HasSubstructMatch(ether_pattern)
    has_carboxyl = mol.HasSubstructMatch(carboxyl_pattern)

    if not (has_hydroxyl or has_epoxy or has_carbonyl or has_ether or has_carboxyl):
        return False, "No hydroxyl, epoxy, carbonyl, ether, or carboxyl groups found"

    return True, "Contains a carotenoid backbone with oxygenated functional groups"