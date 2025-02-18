"""
Classifies: CHEBI:26208 polyunsaturated fatty acid
"""
"""
Classifies: CHEBI:26208 polyunsaturated fatty acid
"""
from rdkit import Chem

def is_polyunsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a polyunsaturated fatty acid based on its SMILES string.
    A polyunsaturated fatty acid is a carboxylic acid with an aliphatic tail containing more than one double bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyunsaturated fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (-COOH)
    carboxylic_acid_pattern = Chem.MolFromSmarts('[CX3](=O)[OX2H1]')
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Count all carbon-carbon double bonds in the molecule
    double_bond_count = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            # Check if both atoms in bond are carbons
            if (bond.GetBeginAtom().GetAtomicNum() == 6 and 
                bond.GetEndAtom().GetAtomicNum() == 6):
                double_bond_count += 1

    if double_bond_count > 1:
        return True, f"Contains {double_bond_count} carbon-carbon double bonds"
    else:
        return False, f"Only {double_bond_count} carbon-carbon double bonds found"