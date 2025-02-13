"""
Classifies: CHEBI:4986 fatty acid methyl ester
"""
from rdkit import Chem

def is_fatty_acid_methyl_ester(smiles: str):
    """
    Determines if a molecule is a fatty acid methyl ester based on its SMILES string.
    A fatty acid methyl ester is characterized by a carbon chain and a methyl ester group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a fatty acid methyl ester, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (None, "Invalid SMILES string")

    # Pattern for methyl ester group: C(=O)OC
    methyl_ester_pattern = Chem.MolFromSmarts("C(=O)OC")
    if not mol.HasSubstructMatch(methyl_ester_pattern):
        return False, "No methyl ester group found"

    # Pattern for a carbon chain: Allow for flexibility in chain fragment length
    # and branching (include terminal carbons at least 6 for smaller ones)
    carbon_chain_pattern = Chem.MolFromSmarts("[CX4]-[CX4]-[CX4]-[CX4]-[CX4]-[CX4]")
    if not mol.HasSubstructMatch(carbon_chain_pattern):
        # Inspect specific conditions where shorter/branched chains qualify
        flexible_chain_pattern = Chem.MolFromSmarts("[CX4]-[CX4]-[CX4]-[C]-[CX4]")
        if not mol.HasSubstructMatch(flexible_chain_pattern):
            return False, "Insufficient/supportive carbon chain length or structure"
    
    return True, "Contains a carbon chain with a methyl ester group"