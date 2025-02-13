"""
Classifies: CHEBI:35819 branched-chain fatty acid
"""
"""
Classifies: CHEBI:38174 branched-chain fatty acid

Any fatty acid in which the parent hydrocarbon chain has one or more alkyl substituents; 
a common component in animal and bacterial lipids. The fatty acyl chain is usually saturated 
and the substituent a methyl group; however, unsaturated BCFAs are found in marine animals, 
and branches other than methyl are found in microbial lipids.
"""
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_branched_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acid based on its SMILES string.

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
    
    # Look for carboxylic acid group (-C(=O)O)
    acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Look for long carbon chain (>= 6 carbons)
    chain_pattern = Chem.MolFromSmarts("[C;H3][C;H2]([C;H2])[C;H2]([C;H2])[C;H2][C;H2]")
    chain_matches = mol.GetSubstructMatches(chain_pattern)
    if not chain_matches:
        return False, "Carbon chain too short for fatty acid"
    
    # Check for alkyl substituents (branches)
    branch_pattern = Chem.MolFromSmarts("[C][C]([C])([C])[C]")
    branch_matches = mol.GetSubstructMatches(branch_pattern)
    if not branch_matches:
        return False, "No alkyl substituents (branches) found"
    
    # Count rotatable bonds, must be >= 2 for branches
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 2:
        return False, "Not enough rotatable bonds for branching"
    
    # Count carbons, must be >= 6
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 6:
        return False, "Too few carbons for fatty acid"
    
    # Count oxygens, must be 2 (carboxylic acid)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count != 2:
        return False, "Incorrect number of oxygens for fatty acid"
    
    return True, "Contains carboxylic acid and long carbon chain with alkyl substituents (branches)"