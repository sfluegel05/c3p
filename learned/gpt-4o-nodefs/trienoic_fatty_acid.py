"""
Classifies: CHEBI:73155 trienoic fatty acid
"""
from rdkit import Chem

def is_trienoic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a trienoic fatty acid based on its SMILES string.
    A trienoic fatty acid typically has a carboxylic acid end group, an aliphatic chain with
    three double bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a trienoic fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Identify double bonds (general case)
    double_bond_pattern = Chem.MolFromSmarts("[CD2]=[CD2]")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    if len(double_bond_matches) < 3:
        return False, f"Found {len(double_bond_matches)} double bonds, need at least 3 for trienoic structure"

    # Check for long carbon chain via carbon count
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 12:
        return False, "Too short carbon chain for a typical trienoic fatty acid"

    # Ensure that double bonds form part of a continuous carbon chain
    # Simplified check: the molecule's carbon backbone should allow for a reasonable chain of carbons
    sp2_carbon_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetHybridization() == Chem.rdchem.HybridizationType.SP2]
    if len(sp2_carbon_indices) < 3:
        return False, "Not enough sp2 hybridized carbons in a chain to indicate trienoic characteristic"

    # Additional flexibility in chain definition
    alkene_chain_pattern = Chem.MolFromSmarts("C=C")
    if not mol.HasSubstructMatch(alkene_chain_pattern):
        return False, "Lacks continuous alkene chain structure typical of trienoic fatty acids"

    return True, "Contains carboxylic acid group and three or more double bonds in a long carbon chain typical of trienoic fatty acid"