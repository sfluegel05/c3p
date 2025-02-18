"""
Classifies: CHEBI:27208 unsaturated fatty acid
"""
from rdkit import Chem

def is_unsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is an unsaturated fatty acid based on its SMILES string.
    An unsaturated fatty acid is a carboxylic acid with an aliphatic chain containing at least one C=C or C#C bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an unsaturated fatty acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (including deprotonated form)
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[O;H1,O-]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Collect atoms in carboxylic acid groups to exclude their double bonds
    carboxylic_atoms = set()
    for match in mol.GetSubstructMatches(carboxylic_acid_pattern):
        # The pattern matches C(=O)O (protonated or deprotonated)
        # Indices: 0=C, 1=O (from =O), 2=O (from -O- or -O-)
        carboxylic_atoms.update(match)

    # Check for any double or triple bonds not in the carboxylic acid group
    has_unsaturation = False
    for bond in mol.GetBonds():
        if bond.GetBondType() in (Chem.BondType.DOUBLE, Chem.BondType.TRIPLE):
            a1 = bond.GetBeginAtomIdx()
            a2 = bond.GetEndAtomIdx()
            # Check if neither atom is part of the carboxylic acid group
            if a1 not in carboxylic_atoms and a2 not in carboxylic_atoms:
                has_unsaturation = True
                break

    if not has_unsaturation:
        return False, "No unsaturation (C=C or C#C) in the aliphatic chain"

    return True, "Carboxylic acid with unsaturation in the aliphatic chain"