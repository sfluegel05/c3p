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

    # Check for disqualifying functional groups (ester, amide, phosphate)
    ester_pattern = Chem.MolFromSmarts('[OD1]-C(=O)-C')
    amide_pattern = Chem.MolFromSmarts('[ND1]-C(=O)')
    phosphate_pattern = Chem.MolFromSmarts('[P](=O)(O)(O)')
    
    if mol.HasSubstructMatch(ester_pattern):
        return False, "Contains ester group"
    if mol.HasSubstructMatch(amide_pattern):
        return False, "Contains amide group"
    if mol.HasSubstructMatch(phosphate_pattern):
        return False, "Contains phosphate group"

    # Count non-aromatic carbon-carbon double bonds
    double_bond_count = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            begin = bond.GetBeginAtom()
            end = bond.GetEndAtom()
            if begin.GetAtomicNum() == 6 and end.GetAtomicNum() == 6:
                # Exclude conjugated systems in non-aliphatic parts
                if not bond.IsInRing() or bond.GetIsConjugated():
                    double_bond_count += 1

    if double_bond_count > 1:
        return True, f"Contains {double_bond_count} carbon-carbon double bonds"
    else:
        return False, f"Only {double_bond_count} carbon-carbon double bonds found"