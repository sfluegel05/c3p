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
    An alkene is an acyclic branched or unbranched hydrocarbon having one carbon-carbon double bond.
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

    # Check that all atoms are carbon or hydrogen
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in (1, 6):
            return False, f"Contains atom other than carbon or hydrogen: {atom.GetSymbol()}"

    # Check for rings involving carbon atoms
    ri = mol.GetRingInfo()
    if ri.NumRings() > 0:
        # Check if any ring contains a carbon atom
        atom_rings = ri.AtomRings()
        for ring in atom_rings:
            if any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring):
                return False, "Molecule contains carbon atoms in rings"

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
            # Exclude triple bonds or other bond types
            return False, f"Contains non-single bond: {bond_type}"

    if num_cc_double_bonds != 1:
        return False, f"Contains {num_cc_double_bonds} carbon-carbon double bonds, requires exactly 1"

    return True, "Molecule is an acyclic hydrocarbon with exactly one carbon-carbon double bond"

__metadata__ = {
    'chemical_class': {
        'name': 'alkene',
        'definition': 'An acyclic branched or unbranched hydrocarbon having one carbon-carbon double bond and the general formula CnH2n. Acyclic branched or unbranched hydrocarbons having more than one double bond are alkadienes, alkatrienes, etc.'
    }
}