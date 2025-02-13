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

    # Check for long aliphatic chain
    chain_pattern = Chem.MolFromSmarts("[CH3][CH2]([CH2])[CH2]([CH2])[CH2]([CH2])[CH2]([CH2])[CH2]([CH2])[CH2]([CH2])[CH2]([CH2])[CH2]([CH2])[CH2]([CH2])[CH2]([CH2])[CH2]([CH2])[CH2]([CH2])[CH2]([CH2])[CH2]([CH2])[CH2]([CH2])[CH2]")
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "No long aliphatic chain found"

    # Check for common double bond positions
    double_bond_patterns = [
        Chem.MolFromSmarts("[CH2]=[CH][CH2][CH2]"),  # cis
        Chem.MolFromSmarts("[CH2]=[CH][CH2][CH3]"),  # trans
        Chem.MolFromSmarts("[CH2]=[CH][CH2][CH2]=[CH][CH2]"),  # cis,cis
        Chem.MolFromSmarts("[CH2]=[CH][CH2][CH3]=[CH][CH2]"),  # cis,trans
        # Add more patterns as needed
    ]
    has_double_bond = any(mol.HasSubstructMatch(pattern) for pattern in double_bond_patterns)

    # Check for common substituents like hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    has_hydroxyl = mol.HasSubstructMatch(hydroxyl_pattern)

    # Count carbon atoms in the longest chain
    longest_chain = rdMolDescriptors.CalcMolLongestChain(mol)
    n_carbon = rdMolDescriptors.CalcLowestChainMultiplier(mol)
    if n_carbon < 13 or n_carbon > 22:
        return False, f"Carbon chain length ({n_carbon}) outside the range of C13 to C22"

    # Check molecular weight - long-chain fatty acids typically >200 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, "Molecular weight too low for long-chain fatty acid"

    # Check number of rotatable bonds - long-chain fatty acids typically have many
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 8:
        return False, "Too few rotatable bonds for long-chain fatty acid"

    reason = "Molecule contains a terminal carboxylic acid group and a long aliphatic chain with a length between C13 and C22"
    if has_double_bond:
        reason += ", with one or more double bonds"
    if has_hydroxyl:
        reason += ", with one or more hydroxyl groups"

    return True, reason