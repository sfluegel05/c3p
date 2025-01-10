"""
Classifies: CHEBI:10615 omega-hydroxy fatty acid
"""
from rdkit import Chem

def is_omega_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an omega-hydroxy fatty acid based on its SMILES string.
    An omega-hydroxy fatty acid has a carboxyl group at position 1 and a hydroxyl group at the last position along a linear chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an omega-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for linear carbon chain with terminal carboxyl and omega hydroxyl group
    # Look for terminal carboxyl group at the end of a chain
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No terminal carboxyl group found"

    # Look for omega hydroxyl group (should be at the last carbon in a linear chain)
    hydroxyl_pattern = Chem.MolFromSmarts("O[CX4]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)

    # Find the linear carbon chain's length by ensuring continuous single bonded carbons
    longest_chain = max([len(path) for path in Chem.rdmolops.GetLongestPath(mol)])
    
    # Ensure the hydroxyl is omega
    for match in hydroxyl_matches:
        hydroxyl_idxs = match
        if len(hydroxyl_idxs) == 2:
            chain_start_idx, hydroxyl_carbon_idx = hydroxyl_idxs
            if longest_chain > 5 and hydroxyl_carbon_idx == longest_chain:
                return True, "Molecule is an omega-hydroxy fatty acid: has terminal carboxyl and omega-hydroxyl groups on a linear chain"

    return False, "Does not conform to omega-hydroxy fatty acid structure"

# Example usage
# smiles = "OCCCCCCCCCCCCCCCCCCCC(O)=O"
# result, reason = is_omega_hydroxy_fatty_acid(smiles)
# print(result, reason)