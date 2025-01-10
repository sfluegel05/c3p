"""
Classifies: CHEBI:39437 tocol
"""
from rdkit import Chem

def is_tocol(smiles: str):
    """
    Determines if a molecule is a tocol based on its SMILES string.
    Tocols have a chromanol core with a hydrocarbon chain attached at position 2
    consisting of three isoprenoid units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tocol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a general chromanol pattern that allows for variations
    chromanol_pattern = Chem.MolFromSmarts("c1cc(O)c2c(c1)CCC2O")
    if not mol.HasSubstructMatch(chromanol_pattern):
        return False, "No chroman-6-ol core found"

    # Define isoprenoid-like pattern allowing for some flexibility in chain
    isoprenoid_pattern = Chem.MolFromSmarts("C(C)(C)CC")
    chains = [match for match in mol.GetSubstructMatches(isoprenoid_pattern)]
    
    # Check if hydrocarbon chain is long enough (at least 15 carbons typical for 3 isoprenoid units)
    hydrocarbon_chains = Chem.MolFromSmarts("[C](C)(C)C")
    long_chain = False
    for match in mol.GetSubstructMatches(hydrocarbon_chains):
        if len(match) >= 15:
            long_chain = True
            break
    
    if not long_chain:
        return False, "Hydrocarbon chain doesn't meet length criteria for 3 isoprenoid units"

    # Ensure correct attachment at position 2 of the chromanol
    core_match = mol.GetSubstructMatch(chromanol_pattern)
    if core_match:
        # Given patterns, verify attachment at correct position
        atom_idx_position_2 = core_match[3]  # index adjustment for chromanol position
        for match in chains:
            if atom_idx_position_2 in match:
                return True, "Contains tocol structure with required chromanol core and hydrocarbon chain"

    return False, "Hydrocarbon chain not correctly attached to chromanol core"