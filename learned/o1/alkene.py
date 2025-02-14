"""
Classifies: CHEBI:32878 alkene
"""
"""
Classifies: alkene
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_alkene(smiles: str):
    """
    Determines if a molecule is an alkene based on its SMILES string.
    An alkene is an acyclic branched or unbranched hydrocarbon having one carbon-carbon double bond
    and the general formula CnH2n. Acyclic branched or unbranched hydrocarbons having more than one
    double bond are alkadienes, alkatrienes, etc.

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

    # Check for acyclic structure
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule is cyclic, alkenes are acyclic"

    # Check that all atoms are carbon or hydrogen
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in (1, 6):
            return False, f"Contains atom other than carbon or hydrogen: {atom.GetSymbol()}"

    # Count number of carbon atoms and hydrogen atoms
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    num_hydrogens = sum(atom.GetTotalNumHs() for atom in mol.GetAtoms())

    # Verify molecular formula matches CnH2n
    if num_hydrogens != 2 * num_carbons:
        return False, f"Molecular formula is C{num_carbons}H{num_hydrogens}, does not match CnH2n"

    # Count number of carbon-carbon double bonds
    num_cc_double_bonds = 0
    for bond in mol.GetBonds():
        begin_atom = bond.GetBeginAtom()
        end_atom = bond.GetEndAtom()
        bond_type = bond.GetBondType()

        if bond_type == Chem.rdchem.BondType.DOUBLE:
            # Check if double bond is between two carbons
            if begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 6:
                num_cc_double_bonds += 1
            else:
                return False, "Contains double bonds involving atoms other than carbon"

        elif bond_type != Chem.rdchem.BondType.SINGLE:
            # Exclude triple bonds or aromatic bonds
            return False, f"Contains non-single/non-double bond: {bond_type}"

    if num_cc_double_bonds != 1:
        return False, f"Contains {num_cc_double_bonds} carbon-carbon double bonds, requires exactly 1"

    return True, "Molecule is an acyclic hydrocarbon with exactly one carbon-carbon double bond and formula CnH2n"

__metadata__ = {
    'chemical_class': {
        'name': 'alkene',
        'definition': 'An acyclic branched or unbranched hydrocarbon having one carbon-carbon double bond and the general formula CnH2n. Acyclic branched or unbranched hydrocarbons having more than one double bond are alkadienes, alkatrienes, etc.'
    }
}