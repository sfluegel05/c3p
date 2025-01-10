"""
Classifies: CHEBI:35819 branched-chain fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_branched_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acid based on its SMILES string.
    A branched-chain fatty acid has a long carbon chain with one or more alkyl substituents and a carboxylic acid group.

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

    # Check for carboxylic acid group (R-C(=O)OH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Check for long alkyl chain (R-CH2-R)
    alkyl_chain_pattern = Chem.MolFromSmarts("[CX3](C)C")
    alkyl_matches = mol.GetSubstructMatches(alkyl_chain_pattern)
    if len(alkyl_matches) < 5:  # Typically long chains in fatty acids
        return False, "Alkyl chain too short"

    # Look for branching (branches must exist for BCFA)
    branch_pattern = Chem.MolFromSmarts("[CX4][CX4]([CX4])C")
    if not mol.HasSubstructMatch(branch_pattern):
        return False, "No branch points found"

    # Estimate molecular weight - emphasizing long chain
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, "Molecular weight too low for typical BCFA"

    return True, "Contains long carbon chain with branch points and a carboxylic acid group"