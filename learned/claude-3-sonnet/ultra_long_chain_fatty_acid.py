"""
Classifies: CHEBI:143004 ultra-long-chain fatty acid
"""
"""
Classifies: CHEBI:51165 ultra-long-chain fatty acid
An ultra-long-chain fatty acid is any very long-chain fatty acid which has a chain length greater than C27.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_ultra_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is an ultra-long-chain fatty acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an ultra-long-chain fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Sanitize molecule
    AllChem.MMFFSanitizeMolecule(mol)
    
    # Check for carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "Does not contain a carboxylic acid group"
    
    # Find the length of the longest carbon chain
    chain_length = 0
    for length in range(28, 100):  # Iterate over chain lengths from 28 to 99
        chain_pattern = Chem.MolFromSmarts(f"[C]~[C]~[C]~[C]~[C]~[C]~[C]~{f'~[C]' * (length - 7)}")
        if mol.HasSubstructMatch(chain_pattern):
            chain_length = length
            break
    
    if chain_length <= 27:
        return False, f"Carbon chain length is only {chain_length}, need greater than 27"
    
    # Check for the presence of only single bonds in the carbon chain
    for bond in mol.GetBonds():
        if bond.GetBondType() != Chem.BondType.SINGLE:
            return False, "Carbon chain contains multiple bonds"
    
    # Check for the absence of additional functional groups
    allowed_atoms = [6, 8]  # Carbon and oxygen
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atoms:
            return False, "Contains additional functional groups besides the carboxylic acid"
    
    return True, f"Contains a single, unbranched carbon chain of length {chain_length} with a carboxylic acid group"