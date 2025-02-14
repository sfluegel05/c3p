"""
Classifies: CHEBI:26208 polyunsaturated fatty acid
"""
from rdkit import Chem

def is_polyunsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a polyunsaturated fatty acid based on its SMILES string.
    A polyunsaturated fatty acid is a fatty acid containing more than one double bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyunsaturated fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a carboxylic acid group (C(=O)O)
    carboxylic_acid = Chem.MolFromSmarts('C(=O)[O;H1,-1]')
    if not mol.HasSubstructMatch(carboxylic_acid):
        return False, "No carboxylic acid group found"

    # Initialize count of aliphatic carbon-carbon double bonds
    double_bond_count = 0

    # Iterate over all bonds in the molecule
    for bond in mol.GetBonds():
        # Check if the bond is a double bond
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            # Get the atoms connected by the bond
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            # Check if both atoms are carbon
            if begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 6:
                # Exclude bonds in rings and aromatic bonds
                if not bond.IsInRing() and not bond.GetIsAromatic():
                    double_bond_count += 1

    # Check if there are more than one double bonds
    if double_bond_count > 1:
        return True, f"Contains {double_bond_count} aliphatic carbon-carbon double bonds"
    else:
        return False, "Does not contain more than one aliphatic carbon-carbon double bond"