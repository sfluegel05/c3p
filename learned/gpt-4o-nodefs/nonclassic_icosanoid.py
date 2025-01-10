"""
Classifies: CHEBI:61703 nonclassic icosanoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nonclassic_icosanoid(smiles: str):
    """
    Determines if a molecule is a nonclassic icosanoid based on its SMILES string.
    Nonclassic icosanoids often include epoxy groups, multiple hydroxyl groups,
    polyunsaturation, and a terminal carboxylic acid. They generally have around 20 carbons.

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
    
    # Flexible epoxy group detection
    epoxy_pattern = Chem.MolFromSmarts("C1OC1")
    epoxy_matches = mol.GetSubstructMatches(epoxy_pattern)
    
    # Polyunsaturations: count overall double bonds
    n_double_bonds = len([bond for bond in mol.GetBonds() if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE])
    
    # Check hydroxyl groups; presence of multiple rather than fixed number
    hydroxyl_pattern = Chem.MolFromSmarts("[CX4][OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    
    # Check for carboxylic acid terminal group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Count the total carbon atoms
    total_carbon_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    # Final evaluation using refined criteria
    if len(epoxy_matches) >= 1 and len(hydroxyl_matches) >= 2 and n_double_bonds >= 3 and total_carbon_atoms >= 18 and total_carbon_atoms <= 22:
        return True, "Matches the criteria for nonclassic icosanoids"

    return False, "Does not satisfy all criteria for a nonclassic icosanoid"