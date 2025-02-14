"""
Classifies: CHEBI:59549 essential fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_essential_fatty_acid(smiles: str):
    """
    Determines if a molecule is an essential fatty acid based on its SMILES string.
    Essential fatty acids are polyunsaturated fatty acids that are required in the diet.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an essential fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carboxylic acid group
    acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxylic acid group found"

    # Count carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 16:
        return False, "Too few carbons for an essential fatty acid"
        
    # Check for double bonds (both cis and trans, as both exist in examples)
    double_bond_pattern = Chem.MolFromSmarts("[CX3]=[CX3]")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    if len(double_bond_matches) < 2:
        return False, f"Found {len(double_bond_matches)} double bonds, need at least 2"


    #check for long aliphatic chain
    aliphatic_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    aliphatic_matches = mol.GetSubstructMatches(aliphatic_chain_pattern)
    if len(aliphatic_matches) < 1:
        return False, f"Missing aliphatic chain, got {len(aliphatic_matches)}"

    # Check the molecule weight, long chains
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250 : #Lower bound based on example molecules
        return False, "Molecular weight is too low for a long chain fatty acid"

    return True, "Contains a carboxylic acid group, long carbon chain, and multiple double bonds."