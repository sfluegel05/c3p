"""
Classifies: CHEBI:32878 alkene
"""
"""
Classifies: alkene
"""
from rdkit import Chem

def is_alkene(smiles: str):
    """
    Determines if a molecule is an alkene based on its SMILES string.
    An alkene is a hydrocarbon (acyclic or cyclic, branched or unbranched) having one carbon-carbon double bond.
    Acyclic branched or unbranched hydrocarbons having more than one double bond are alkadienes, alkatrienes, etc.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkene, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count number of carbon atoms
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons == 0:
        return False, "No carbon atoms found"

    # Count number of carbon-carbon double bonds
    num_cc_double_bonds = 0
    for bond in mol.GetBonds():
        bond_type = bond.GetBondType()
        if bond_type == Chem.rdchem.BondType.DOUBLE:
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            if begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 6:
                num_cc_double_bonds += 1
        elif bond_type == Chem.rdchem.BondType.TRIPLE:
            return False, "Contains triple bonds, not an alkene"
        elif bond_type != Chem.rdchem.BondType.SINGLE:
            return False, f"Contains non-single/non-double bond: {bond_type}"

    if num_cc_double_bonds != 1:
        return False, f"Contains {num_cc_double_bonds} carbon-carbon double bonds, requires exactly 1"

    return True, "Molecule has exactly one carbon-carbon double bond"

__metadata__ = {
    'chemical_class': {
        'name': 'alkene',
        'definition': 'An acyclic branched or unbranched hydrocarbon having one carbon-carbon double bond and the general formula CnH2n. Acyclic branched or unbranched hydrocarbons having more than one double bond are alkadienes, alkatrienes, etc.'
    }
}