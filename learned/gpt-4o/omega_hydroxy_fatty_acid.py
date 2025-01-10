"""
Classifies: CHEBI:10615 omega-hydroxy fatty acid
"""
from rdkit import Chem

def is_omega_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an omega-hydroxy fatty acid based on its SMILES string.
    An omega-hydroxy fatty acid has a carboxyl group at position 1 and a hydroxyl at the
    last carbon atom of a straight-chain fatty acid.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is an omega-hydroxy fatty acid, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify main carbon chain candidates (long carbon sequences)
    carbon_chain_pattern = Chem.MolFromSmarts("C[C,C](C)[C,C](C)")
    chain_candidates = mol.GetSubstructMatches(carbon_chain_pattern)

    if not chain_candidates:
        return False, "No suitable carbon chain detected"

    # Check each candidate for terminal carboxyl and hydroxyl
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    hydroxyl_pattern = Chem.MolFromSmarts("[CX4;H2,H3][OH]")

    for chain in chain_candidates:
        start_atom = mol.GetAtomWithIdx(chain[0])
        end_atom = mol.GetAtomWithIdx(chain[-1])

        # Check for carboxyl group at start
        carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
        if not any(start_atom.GetIdx() in match for match in carboxyl_matches):
            continue

        # Check for hydroxyl group at end
        hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
        if not any(end_atom.GetIdx() == match[1] for match in hydroxyl_matches):
            continue

        # If both checks pass
        return True, "Contains carboxyl group at position 1 and omega-hydroxyl group at end of the main chain as per omega-hydroxy fatty acid definition."

    return False, "No omega-hydroxy fatty acid structure detected"