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

    # Find the carbon atom in the carboxylic acid group
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if not carboxylic_acid_matches:
        return False, "No carboxylic acid group found"
    
    # Get the carbon atom in the carboxylic acid group
    carboxylic_carbon = carboxylic_acid_matches[0][0]

    # Check for a hydroxyl group (-OH) at the 2-position (alpha position)
    alpha_carbon = None
    for neighbor in mol.GetAtomWithIdx(carboxylic_carbon).GetNeighbors():
        if neighbor.GetAtomicNum() == 6:  # Carbon atom
            alpha_carbon = neighbor.GetIdx()
            break

    if alpha_carbon is None:
        return False, "No alpha carbon found"

    # Check if the alpha carbon has a hydroxyl group
    alpha_carbon_atom = mol.GetAtomWithIdx(alpha_carbon)
    has_hydroxyl = False
    for neighbor in alpha_carbon_atom.GetNeighbors():
        if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() == 1:  # Hydroxyl group
            has_hydroxyl = True
            break

    if not has_hydroxyl:
        return False, "No hydroxyl group at the 2-position"

    # Check for a sufficiently long carbon chain (fatty acid)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 6:
        return False, "Carbon chain too short to be a fatty acid"

    return True, "Contains a carboxylic acid group, a hydroxyl group at the 2-position, and a sufficiently long carbon chain"