from rdkit import Chem
from rdkit.Chem import AllChem

def is_sulfonium_betaine(smiles: str):
    """
    Determines if a molecule is a sulfonium betaine - a neutral molecule with charge-separated form
    containing a sulfonium atom with no hydrogens that is not adjacent to the anionic atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sulfonium betaine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Check for presence of sulfonium cation (S+)
    sulfonium_atoms = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'S' and atom.GetFormalCharge() == 1:
            sulfonium_atoms.append(atom)

    if not sulfonium_atoms:
        return False, "No sulfonium cation (S+) found"

    # Check if any sulfonium has hydrogens
    for s_atom in sulfonium_atoms:
        if s_atom.GetTotalNumHs() > 0:
            return False, "Sulfonium has hydrogen atoms"

    # Find anionic atoms (atoms with negative charge)
    anionic_atoms = []
    for atom in mol.GetAtoms():
        if atom.GetFormalCharge() == -1:
            anionic_atoms.append(atom)

    if not anionic_atoms:
        return False, "No anionic atoms found"

    # Check if molecule is overall neutral
    total_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    if total_charge != 0:
        return False, "Molecule is not neutral overall"

    # Check if any sulfonium is adjacent to anionic atom
    for s_atom in sulfonium_atoms:
        s_neighbors = [n.GetIdx() for n in s_atom.GetNeighbors()]
        for a_atom in anionic_atoms:
            if a_atom.GetIdx() in s_neighbors:
                return False, "Sulfonium is adjacent to anionic atom"

    return True, "Molecule is a sulfonium betaine"
# Pr=1.0
# Recall=1.0