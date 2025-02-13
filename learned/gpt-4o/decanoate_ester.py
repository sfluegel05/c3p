"""
Classifies: CHEBI:87658 decanoate ester
"""
from rdkit import Chem

def is_decanoate_ester(smiles: str):
    """
    Determines if a molecule is a decanoate ester based on its SMILES string.
    A decanoate ester is formed by esterification of decanoic acid with an alcohol or phenol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a decanoate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Decanoic acid chain pattern: 10 carbon atoms in a row
    decanoic_acid_chain_pattern = Chem.MolFromSmarts("CCCCCCCCCC")
    if not mol.HasSubstructMatch(decanoic_acid_chain_pattern):
        return False, "No decanoic acid chain found"

    # Ester linkage pattern
    ester_linkage_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(ester_linkage_pattern):
        return False, "No ester linkage found"

    # Check if an ester group is linked to a decanoic acid chain
    ester_matches = mol.GetSubstructMatches(ester_linkage_pattern)
    for match in ester_matches:
        carbon_index = match[0]  # The carbon in the ester group is the first in the match
        atom = mol.GetAtomWithIdx(carbon_index)
        if any(atom.GetNeighbors()[i].HasSubstructMatch(decanoic_acid_chain_pattern) for i in range(len(atom.GetNeighbors()))):
            return True, "Contains decanoate ester linkage"
    
    return False, "No ester linkage directly connected to a decanoic acid chain"