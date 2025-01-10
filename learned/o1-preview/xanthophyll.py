"""
Classifies: CHEBI:27325 xanthophyll
"""
"""
Classifies: xanthophyll
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_xanthophyll(smiles: str):
    """
    Determines if a molecule is a xanthophyll based on its SMILES string.
    A xanthophyll is an oxygenated carotenoid, consisting of oxygenated carotenes.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a xanthophyll, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of oxygen atoms
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count == 0:
        return False, "No oxygen atoms found; not an oxygenated carotene"

    # Check for a long conjugated chain (number of conjugated double bonds)
    # Approximate by counting the number of conjugated double bonds between carbons
    num_conj_double_bonds = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and bond.GetIsConjugated():
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            if begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 6:
                num_conj_double_bonds += 1

    if num_conj_double_bonds < 9:
        return False, f"Insufficient conjugated double bonds ({num_conj_double_bonds}); not a carotenoid"

    # Check that the molecule has sufficient carbon atoms (approximate carotenoid size)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 30:
        return False, f"Too few carbon atoms ({c_count}); not a carotenoid"

    # Check molecular weight - carotenoids typically >500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, f"Molecular weight too low ({mol_wt:.2f} Da); not a carotenoid"

    return True, "Molecule is an oxygenated carotenoid (xanthophyll)"

__metadata__ = {   
    'chemical_class': {   
        'id': 'CHEBI:27530',
        'name': 'xanthophyll',
        'definition': 'A subclass of carotenoids consisting of the oxygenated carotenes.',
    },
}