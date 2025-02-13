"""
Classifies: CHEBI:15904 long-chain fatty acid
"""
"""
Classifies: ChEBI:38898 Long-chain fatty acid
A fatty acid with a chain length ranging from C13 to C22.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for terminal carboxylic acid group
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No terminal carboxylic acid group found"

    # Check for long unbranched aliphatic chain
    chain_pattern = Chem.MolFromSmarts("[CH3][CH2]([CH2])[CH2]([CH2])[CH2]([CH2])[CH2]([CH2])[CH2]([CH2])[CH2]([CH2])[CH2]([CH2])[CH2]([CH2])[CH2]([CH2])[CH2]([CH2])[CH2]([CH2])[CH2]([CH2])[CH2]([CH2])[CH2]([CH2])[CH2]([CH2])[CH2]([CH2])[CH2]([CH2])[CH2]([CH2])[CH2]([CH2])[CH2]([CH2])[CH2](C(=O)O)[CH2][CH3]")
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "No long unbranched aliphatic chain found"

    # Check for double bonds at specific positions (optional)
    # ...

    # Count carbon atoms in the longest chain
    longest_chain = rdMolDescriptors.CalcMolLongestChain(mol)
    n_carbon = rdMolDescriptors.CalcLowestChainMultiplier(mol)
    if n_carbon < 13 or n_carbon > 22:
        return False, f"Carbon chain length ({n_carbon}) outside the range of C13 to C22"

    return True, "Molecule contains a terminal carboxylic acid group and a long unbranched aliphatic chain with a length between C13 and C22"