"""
Classifies: CHEBI:61910 very long-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_very_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a very long-chain fatty acyl-CoA based on its SMILES string.
    A very long-chain fatty acyl-CoA has a fatty acyl group with a chain length greater than C22.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a very long-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for coenzyme A substructure
    # This SMARTS pattern might need adjustment based on the examples provided
    coa_pattern = Chem.MolFromSmarts("C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)C")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No Coenzyme A backbone found"

    # Identify the longest carbon chain
    # Collect all carbon atoms and determine the longest path based on position
    atom_count_map = {}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon atom
            neighbors = [x.GetIdx() for x in atom.GetNeighbors() if x.GetAtomicNum() == 6]
            atom_count_map[atom.GetIdx()] = len(neighbors)

    # Calculate longest path
    if atom_count_map:
        longest_chain_length = max(atom_count_map.values())
        if longest_chain_length <= 22:
            return False, f"Longest carbon chain length is {longest_chain_length}, not greater than C22"
    else:
        return False, "No carbon chain found"

    return True, f"Contains CoA backbone and fatty acyl chain length is {longest_chain_length}, which is greater than C22"