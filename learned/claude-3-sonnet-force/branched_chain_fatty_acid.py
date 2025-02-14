"""
Classifies: CHEBI:35819 branched-chain fatty acid
"""
"""
Classifies: CHEBI:32857 branched-chain fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_branched_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acid based on its SMILES string.
    A branched-chain fatty acid is a fatty acid with one or more alkyl substituents on the parent hydrocarbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a branched-chain fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carboxylic acid group
    acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Check for alkyl branches
    alkyl_pattern = Chem.MolFromSmarts("[CX4H3]")
    branch_matches = mol.GetSubstructMatches(alkyl_pattern)
    if len(branch_matches) < 1:
        return False, "No alkyl branches found"
    
    # Check for long carbon chain (at least 4 carbons)
    longest_chain = AllChem.GetLongestPathWithMultipleTraversals(mol)
    if len(longest_chain) < 4:
        return False, "Carbon chain too short for fatty acid"
    
    # Count hydrogens to ensure saturated or unsaturated (no cycles/aromatics)
    h_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 1)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if h_count != (2*c_count + 2):
        return False, "Contains cycles or aromatic rings"
    
    # Check for halogen substituents
    halogen_pattern = Chem.MolFromSmarts("[F,Cl,Br,I]")
    if mol.HasSubstructMatch(halogen_pattern):
        return False, "Contains halogen substituents"
    
    # All checks passed
    return True, "Contains a carboxylic acid group and alkyl branches on a long hydrocarbon chain"