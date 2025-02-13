"""
Classifies: CHEBI:25413 monounsaturated fatty acid
"""
"""
Classifies: CHEBI:36975 monounsaturated fatty acid
Any fatty acid with one double or triple bond in the fatty acid chain and singly bonded carbon atoms in the rest of the chain. MUFAs have positive effects on the cardiovascular system, and in diabetes treatment.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_monounsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monounsaturated fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Check for carbon chain with at least 4 carbons
    chain_pattern = Chem.MolFromSmarts("[C;H3]-[C;H2]~[C;H2]~[C;H2]")
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "Carbon chain too short"
    
    # Check for exactly one double or triple bond, ignoring oxygen-containing groups
    double_bond_pattern = Chem.MolFromSmarts("[C;H2]=[C;H2]")
    triple_bond_pattern = Chem.MolFromSmarts("[C;H2]#[C;H2]")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern, useOpReversedChains=True, useQueryRootDescriptors=True, ignoreOxidePermuteDerivatives=True)
    triple_bond_matches = mol.GetSubstructMatches(triple_bond_pattern, useOpReversedChains=True, useQueryRootDescriptors=True, ignoreOxidePermuteDerivatives=True)
    
    if len(double_bond_matches) + len(triple_bond_matches) != 1:
        return False, "Found multiple or no double/triple bonds"
    
    # Check for reasonable molecular weight and atom counts
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 100 or mol_wt > 700:
        return False, "Molecular weight outside typical range for fatty acids"
    
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 4 or o_count > 4:
        return False, "Unusual atom counts for a fatty acid"
    
    # Exclude deuterated compounds (optional)
    if any(atom.GetIsotope() > 0 for atom in mol.GetAtoms()):
        return False, "Deuterated compound detected"
    
    return True, "Contains exactly one double or triple bond in a carbon chain with a carboxylic acid group"