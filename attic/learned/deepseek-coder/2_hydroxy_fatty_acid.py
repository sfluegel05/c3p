"""
Classifies: CHEBI:10283 2-hydroxy fatty acid
"""
"""
Classifies: 2-hydroxy fatty acid
Definition: Any fatty acid with a hydroxy functional group in the alpha- or 2-position.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_2_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxy fatty acid based on its SMILES string.
    A 2-hydroxy fatty acid has a carboxylic acid group and a hydroxyl group at the 2-position,
    along with a sufficiently long carbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (-C(=O)OH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Find the carboxylic acid group
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if not carboxylic_acid_matches:
        return False, "No carboxylic acid group found"

    # Iterate over all carboxylic acid groups to find one with a hydroxyl at the 2-position
    for match in carboxylic_acid_matches:
        carboxylic_acid_atom = match[0]  # The carbon atom of the carboxylic acid group

        # Get the atom connected to the carboxylic acid carbon (the alpha carbon)
        alpha_carbon = None
        for neighbor in mol.GetAtomWithIdx(carboxylic_acid_atom).GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # Carbon atom
                alpha_carbon = neighbor.GetIdx()
                break

        if alpha_carbon is None:
            continue  # No alpha carbon found, skip to the next match

        # Check if the alpha carbon has a hydroxyl group
        for neighbor in mol.GetAtomWithIdx(alpha_carbon).GetNeighbors():
            if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() == 1:  # Hydroxyl group
                # Check for a sufficiently long carbon chain (fatty acid)
                carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
                if carbon_count >= 4:
                    return True, "Contains a carboxylic acid group, a hydroxyl group at the 2-position, and a sufficiently long carbon chain"

    return False, "No hydroxyl group at the 2-position relative to the carboxylic acid group"