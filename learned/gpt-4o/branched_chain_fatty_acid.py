"""
Classifies: CHEBI:35819 branched-chain fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_branched_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is a branched-chain fatty acid, False otherwise.
        str: Reason for classification.
    """
    
    # Parse the SMILES string into an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Check for a terminal carboxylic acid group.
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No terminal carboxylic acid group found"

    # Check the chain length (number of carbon atoms in the longest alkyl chain).
    num_carbons = sum(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms())
    if num_carbons < 7:
        return False, "Too few carbon atoms for a fatty acid (typical count is 7-28)"
    
    # Define SMARTS patterns for common branching found in fatty acids.
    branching_patterns = [
        Chem.MolFromSmarts("[CX4](C)(C)C"),  # Tertiary carbons
        Chem.MolFromSmarts("[CX4;!R](C)(C)C"),  # Avoid rings, branching off-ring
        Chem.MolFromSmarts("[C](C)(C)CO"),  # Hydroxyl group branching
        Chem.MolFromSmarts("C(C)(C)C=C"),  # Unsaturated branch
    ]
    
    # Check if there is a branching pattern match indicative of a branched-chain.
    for pattern in branching_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains terminal carboxylic acid group and appropriate branching pattern"

    return False, "No suitable alkyl branching pattern found"