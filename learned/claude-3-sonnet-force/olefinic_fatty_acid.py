"""
Classifies: CHEBI:53339 olefinic fatty acid
"""
"""
Classifies: CHEBI:36346 olefinic fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_olefinic_fatty_acid(smiles: str):
    """
    Determines if a molecule is an olefinic fatty acid based on its SMILES string.
    An olefinic fatty acid is defined as any fatty acid containing at least one C=C double bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an olefinic fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"

    # Check for aliphatic chain
    aliphatic_chain_pattern = Chem.MolFromSmarts("[CH3,CH2,CH]~[CH3,CH2,CH]")
    if not mol.HasSubstructMatch(aliphatic_chain_pattern):
        return False, "No aliphatic chain found"

    # Check for double bonds
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    if not double_bond_matches:
        return False, "No double bonds found"

    # Additional checks (optional)
    # Check for common functional groups (hydroxy, epoxy, hydroperoxy, etc.)
    hydroxy_pattern = Chem.MolFromSmarts("O")
    epoxy_pattern = Chem.MolFromSmarts("C1OC1")
    hydroperoxy_pattern = Chem.MolFromSmarts("OO")
    
    # Check molecular weight or carbon count range (optional)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Typical ranges for fatty acids (adjust as needed)
    if mol_wt < 150 or mol_wt > 500:
        return False, "Molecular weight outside typical range for fatty acids"
    if c_count < 10 or c_count > 30:
        return False, "Carbon count outside typical range for fatty acids"

    return True, "Contains carboxylic acid group, aliphatic chain, and at least one C=C double bond"