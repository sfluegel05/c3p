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

    # Look for terminal carboxyl group
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=O)[OX1H0-,OX2H1]")
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    
    if not carboxyl_matches:
        return False, "No terminal carboxyl group detected"

    # Check for hydroxyl group not part of carboxyl, potentially omega position
    hydroxy_pattern = Chem.MolFromSmarts("[CX4H2,CX4H1][OX2H1,OX1H0-]")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)

    if not hydroxy_matches:
        return False, "No terminal hydroxyl group detected"

    # Determine the hydroxyl and carboxyl groups proximity to potential main chain
    for hydroxyl_match in hydroxy_matches:
        hydroxyl_idx = hydroxyl_match[1]  # The -OH group atom index
        for carboxyl_match in carboxyl_matches:
            carboxyl_idx = carboxyl_match[0]  # The C=O group atom index of carboxyl

            # Ensure -OH is at the end of a carbon chain and separate from carboxyl
            if mol.GetBondBetweenAtoms(carboxyl_idx, hydroxyl_idx) is not None:
                continue  # Skip if directly bonded (should not normally happen)

            # Check distance between carboxyl and hydroxyl to assume a "main chain"
            shortest_path = Chem.rdmolops.GetShortestPath(mol, carboxyl_idx, hydroxyl_idx)

            if len(shortest_path) > 4:  # decent length for a fatty acid
                return True, "Detected omega-hydroxy group distant from carboxyl group as per definition."

    return False, "No omega-hydroxy fatty acid structure detected"