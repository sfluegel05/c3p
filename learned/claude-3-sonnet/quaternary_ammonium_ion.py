"""
Classifies: CHEBI:35267 quaternary ammonium ion
"""
"""
Classifies: CHEBI:35485 quaternary ammonium ion
A derivative of ammonium, NH4(+), in which all four of the hydrogens bonded to nitrogen
have been replaced with univalent (usually organyl) groups.
"""
from rdkit import Chem
from rdkit.Chem import rdqueries

def is_quaternary_ammonium_ion(smiles: str):
    """
    Determines if a molecule is a quaternary ammonium ion based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a quaternary ammonium ion, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find quaternary nitrogen atoms
    quat_n_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7
                    and atom.GetFormalCharge() == 1
                    and sum(1 for bond in atom.GetBonds() if bond.GetBondType() == Chem.BondType.SINGLE) == 4]

    if not quat_n_atoms:
        return False, "No quaternary nitrogen atoms found"

    # Check if all substituents are univalent (organyl) groups
    for n_atom in quat_n_atoms:
        for bond in n_atom.GetBonds():
            neighbor = bond.GetOtherAtom(n_atom)
            if neighbor.GetAtomicNum() not in [1, 6, 7, 8, 9, 15, 16, 17, 35]:
                return False, f"Substituent on quaternary nitrogen is not a univalent group: {Chem.GetAtomMoniker(neighbor)}"

    return True, "Contains quaternary nitrogen atom with four univalent (organyl) substituents"