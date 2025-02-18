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

    # Check for exactly one carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts('[CX3](=O)[OX2H1]')
    carboxylic_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if len(carboxylic_matches) != 1:
        return False, f"Found {len(carboxylic_matches)} carboxylic acid groups (need exactly 1)"

    # Check for aromatic rings
    if any(atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return False, "Contains aromatic rings"

    # Check for ester groups (R-O-C(=O)-R') excluding the carboxylic acid
    ester_pattern = Chem.MolFromSmarts('[OX2]-C(=O)-[OX2]')
    if mol.HasSubstructMatch(ester_pattern):
        return False, "Contains ester group"

    # Count all aliphatic carbon-carbon double bonds
    double_bond_count = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            # Check if both atoms are carbons and not in aromatic system
            if (bond.GetBeginAtom().GetAtomicNum() == 6 and 
                bond.GetEndAtom().GetAtomicNum() == 6 and
                not bond.IsInRing() and
                not bond.GetIsAromatic()):
                double_bond_count += 1

    if double_bond_count > 1:
        return True, f"Contains {double_bond_count} aliphatic carbon-carbon double bonds"
    else:
        return False, f"Only {double_bond_count} aliphatic carbon-carbon double bonds found"