"""
Classifies: CHEBI:143004 ultra-long-chain fatty acid
"""
"""
Classifies: ultra-long-chain fatty acid (chain length > C27)
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_ultra_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is an ultra-long-chain fatty acid (chain length > C27).
    Must have a carboxylic acid group and a continuous carbon chain >27 atoms.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Verify carboxylic acid group exists
    acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxylic acid group"

    # Find longest carbon chain using RDKit's method
    main_chain = rdMolDescriptors.FindLongestChain(mol)
    chain_length = len(main_chain)

    # Check if chain includes the carboxylic acid's carbonyl carbon
    acid_atoms = {a for match in mol.GetSubstructMatches(acid_pattern) for a in match}
    chain_contains_acid = any(a in acid_atoms for a in main_chain)

    if chain_contains_acid and chain_length > 27:
        return True, f"Main chain length {chain_length} > C27"
    
    return False, f"Main chain length {chain_length} ≤ C27" if chain_contains_acid else "No valid chain connected to acid group"