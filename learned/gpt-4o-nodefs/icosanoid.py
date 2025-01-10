"""
Classifies: CHEBI:23899 icosanoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_icosanoid(smiles: str):
    """
    Determines if a molecule is an icosanoid based on its SMILES string.
    Icosanoids are derived from 20-carbon arachidonic acid and include structures
    with multiple double bonds and functional groups like hydroxyls and carboxylic acids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an icosanoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for long carbon chain (at least 20 carbons)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:
        return False, "Does not have long carbon chain (at least 20 carbons)"
    
    # Check for multiple double bonds
    double_bond_count = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE)
    if double_bond_count < 3:
        return False, "Insufficient double bonds, requires at least 3"
    
    # Check for functional groups: carboxylic acid or hydroxyl
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    hydroxyl_pattern = Chem.MolFromSmarts("CO")
    
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "Missing carboxylic acid group"
    
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "Missing hydroxyl group"
    
    return True, "Contains long carbon chain with multiple double bonds and common functional groups of icosanoids"