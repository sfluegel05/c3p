"""
Classifies: CHEBI:28963 amino sugar
"""
from rdkit import Chem
from rdkit.Chem import rdqueries

def is_amino_sugar(smiles: str):
    """
    Determines if a molecule is an amino sugar based on its SMILES string.
    An amino sugar is characterized by one or more amino groups replacing hydroxyl groups on a sugar molecule.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amino sugar, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS for amino groups
    amino_group_pattern = Chem.MolFromSmarts("[NX3;H2,H1,H0]")  # Any amine (primary, secondary, tertiary)

    # Check for amino groups
    if not mol.HasSubstructMatch(amino_group_pattern):
        return False, "No amino groups found"

    # Define more inclusive SMARTS pattern for sugar rings like pyranose/furanose
    sugar_ring_pattern = Chem.MolFromSmarts("C1(CO)C(O)C(CO)C(O)OC1 |R|")

    # Check for sugar-like structure
    if not mol.HasSubstructMatch(sugar_ring_pattern):
        return False, "No sugar-like structure found"

    # Make sure amino groups replace hydroxyl groups in sugar structure
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'N':
            # Inspect bonded atoms to check if amino is replacing OH
            nitrogen_neighbors = [n.GetAtomicNum() for n in atom.GetNeighbors()]
            if 6 in nitrogen_neighbors:  # Ensure nitrogen is bonded to a carbon
                carbon_neighbors = sum(1 for a in nitrogen_neighbors if a == 6)
                # An amino sugar should be part of a sugar ring where OH groups can be replaced
                if carbon_neighbors >= 1:
                    # Ensure it is indeed replacing an OH group
                    neighboring_carbon = next(n for n in atom.GetNeighbors() if n.GetSymbol() == 'C')
                    neighbors = [a.GetAtomicNum() for a in neighboring_carbon.GetNeighbors()]
                    if 8 not in neighbors:  # Absence of hydroxyl (oxygen) neighbor confirms replacement
                        return True, "Contains a sugar backbone with amino groups substituting hydroxyl groups."

    return False, "No correct hydroxyl replacement was found by an amino group."