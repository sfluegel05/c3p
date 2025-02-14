"""
Classifies: CHEBI:35366 fatty acid
"""
"""
Classifies: CHEBI:36976 fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acid(smiles: str):
    """
    Determines if a molecule is a fatty acid based on its SMILES string.
    A fatty acid is an aliphatic monocarboxylic acid with an unbranched chain of 4-28 carbons.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for carboxylic acid group (-C(=O)O)
    acid_pattern = Chem.MolFromSmarts("C(=O)[O;H,-]")
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Look for long unbranched carbon chain
    chain_pattern = Chem.MolFromSmarts("[C;H3][C;H2]~[C;H2]~[C;H2][C;H3]")
    chain_matches = mol.GetSubstructMatches(chain_pattern)
    
    if not chain_matches:
        return False, "No long unbranched carbon chain found"
    
    # Count carbons in longest chain
    longest_chain = max([len(AllChem.FindAllPathsOfLengthN(mol, 4, path[0], path[-1])) for path in chain_matches], key=len)
    n_carbons = len(longest_chain) + 1  # +1 to account for carboxylic acid carbon
    
    if n_carbons < 4 or n_carbons > 28:
        return False, f"Chain length of {n_carbons} carbons is outside typical fatty acid range (4-28)"
    
    # Check for branching
    if any(atom.GetDegree() > 3 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6):
        return False, "Carbon chain is branched"
    
    return True, f"Contains unbranched aliphatic carbon chain of length {n_carbons} with carboxylic acid group"