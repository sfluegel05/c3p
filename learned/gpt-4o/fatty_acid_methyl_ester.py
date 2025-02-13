"""
Classifies: CHEBI:4986 fatty acid methyl ester
"""
from rdkit import Chem

def is_fatty_acid_methyl_ester(smiles: str):
    """
    Determines if a molecule is a fatty acid methyl ester based on its SMILES string.
    A fatty acid methyl ester is characterized by a long carbon chain and a methyl ester functional group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acid methyl ester, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Pattern for methyl ester group: C(=O)OC
    methyl_ester_pattern = Chem.MolFromSmarts("C(=O)OC")
    if not mol.HasSubstructMatch(methyl_ester_pattern):
        return False, "No methyl ester group found"

    # Pattern for a long hydrocarbon chain, minimum of 8 contiguous carbons is a common cutoff for fatty acids
    long_chain_pattern = Chem.MolFromSmarts("[C;H2,H3][C;H2,H3][C;H2,H3][C;H2,H3][C;H2,H3][C;H2,H3][C;H2,H3][C;H2,H3]") 
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "Insufficient carbon chain length for fatty acid"

    # Check for oxygen count; must be at least 2 (for ester and possibly other functional groups)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 2:
        return False, "Molecule lacks sufficient oxygen atoms to be a methyl ester"

    return True, "Contains a long carbon chain with a methyl ester group"