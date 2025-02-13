"""
Classifies: CHEBI:37142 organoiodine compound
"""
"""
Classifies: CHEBI:27380 organoiodine compound

An organoiodine compound is a compound containing at least one carbon-iodine bond.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_organoiodine_compound(smiles: str):
    """
    Determines if a molecule is an organoiodine compound based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organoiodine compound, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if molecule contains iodine
    if not any(atom.GetAtomicNum() == 53 for atom in mol.GetAtoms()):
        return False, "No iodine atoms present"

    # Check for carbon-iodine bonds
    c_i_bonds = any(bond.GetBeginAtom().GetAtomicNum() == 6 and
                    bond.GetEndAtom().GetAtomicNum() == 53 for bond in mol.GetBonds())
    if not c_i_bonds:
        return False, "No carbon-iodine bonds found"

    # Check for organic moiety
    organic_moiety = True
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 53:  # Iodine atom
            if not any(neighbor.GetAtomicNum() == 6 for neighbor in atom.GetNeighbors()):
                organic_moiety = False
                break

    if not organic_moiety:
        return False, "Iodine atom not part of an organic moiety"

    # Handle charged/ionic species
    mol_ions = [m for m in AllChem.IonizeMolecule(mol, ionize_everything=True)]
    for ion_mol in mol_ions:
        if is_organoiodine_compound(Chem.MolToSmiles(ion_mol))[0]:
            return True, "Ionized form of an organoiodine compound"

    # Handle tautomers
    tautomers = Chem.MolToSmarts(Chem.MolFromSmiles(smiles), isomericSmiles=True)
    for tautomer in tautomers:
        if is_organoiodine_compound(Chem.MolToSmiles(tautomer))[0]:
            return True, "Tautomeric form of an organoiodine compound"

    return True, "Contains at least one carbon-iodine bond in an organic moiety"