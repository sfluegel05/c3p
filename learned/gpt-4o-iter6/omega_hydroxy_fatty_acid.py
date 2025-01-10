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

    # Look for a terminal carboxyl group
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[O;H1]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No terminal carboxyl group found"

    # Look for a terminal hydroxyl group
    # Ensure it is connected to a carbon chain, not directly after the carboxyl
    hydroxyl_pattern = Chem.MolFromSmarts("[CX4][OH]")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No omega hydroxyl group found"

    # Check the carbon chain length is at least 8
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 8:
        return False, f"Carbon chain too short: {c_count} carbons found"

    return True, "Molecule is an omega-hydroxy fatty acid: has terminal carboxyl and omega-hydroxyl groups on a linear chain"

# Example usage
# smiles = "OCCCCCCCCCCCCCCCCCCCC(O)=O"
# result, reason = is_omega_hydroxy_fatty_acid(smiles)
# print(result, reason)