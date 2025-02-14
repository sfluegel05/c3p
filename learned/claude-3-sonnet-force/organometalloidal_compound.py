"""
Classifies: CHEBI:143084 organometalloidal compound
"""
"""
Classifies: CHEBI:35498 organometalloidal compound

A compound having bonds between one or more metalloid atoms and one or more carbon atoms of an organyl group.
"""

from rdkit import Chem

# List of metalloid atomic numbers
METALLOID_ATOMIC_NUMBERS = [5, 14, 32, 33, 51, 52, 83, 84]  # B, Si, Ge, As, Sb, Te, Bi, Po

def is_organometalloidal_compound(smiles: str):
    """
    Determines if a molecule is an organometalloidal compound based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organometalloidal compound, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for metalloid atoms
    has_metalloid = any(atom.GetAtomicNum() in METALLOID_ATOMIC_NUMBERS for atom in mol.GetAtoms())
    if not has_metalloid:
        return False, "No metalloid atoms present"

    # Check for bonds between metalloid and carbon atoms
    for bond in mol.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        if (atom1.GetAtomicNum() in METALLOID_ATOMIC_NUMBERS and atom2.GetAtomicNum() == 6) or \
           (atom2.GetAtomicNum() in METALLOID_ATOMIC_NUMBERS and atom1.GetAtomicNum() == 6):
            return True, "Contains bond(s) between metalloid and carbon atoms"

    return False, "No bond found between metalloid and carbon atoms"