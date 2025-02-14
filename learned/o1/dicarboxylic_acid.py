"""
Classifies: CHEBI:35692 dicarboxylic acid
"""
"""
Classifies: CHEBI:35692 dicarboxylic acid
"""

from rdkit import Chem

def is_dicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a dicarboxylic acid based on its SMILES string.
    A dicarboxylic acid is any carboxylic acid containing two free carboxy groups (-COOH),
    not involved in ester or amide bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dicarboxylic acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define free carboxylic acid group pattern (-C(=O)OH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O;H1]")

    if carboxylic_acid_pattern is None:
        return False, "Failed to create carboxylic acid group pattern"

    # Find all free carboxylic acid groups in the molecule
    carboxyl_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    num_carboxyl_groups = len(carboxyl_matches)

    # Check if there are exactly two free carboxylic acid groups
    if num_carboxyl_groups == 2:
        # Verify that the carboxyl carbons are not involved in amide or ester bonds
        for match in carboxyl_matches:
            carboxyl_carbon_idx = match[0]
            carboxyl_carbon = mol.GetAtomWithIdx(carboxyl_carbon_idx)
            # Get neighbors of carboxyl carbon excluding the double-bonded oxygen and hydroxyl oxygen
            neighbors = [a for a in carboxyl_carbon.GetNeighbors() if a.GetAtomicNum() != 8]
            # The carboxyl carbon should have only one neighbor that is not oxygen (the attached carbon)
            if len(neighbors) != 1:
                return False, "Carboxyl group is part of an ester or amide"
            # Check if the neighbor is not a nitrogen (amide) or oxygen (ester)
            neighbor_atom_num = neighbors[0].GetAtomicNum()
            if neighbor_atom_num == 7 or neighbor_atom_num == 8:
                return False, "Carboxyl group is part of an ester or amide"

        return True, "Molecule contains exactly two free carboxyl groups"
    else:
        return False, f"Found {num_carboxyl_groups} free carboxyl group(s), expected exactly 2"