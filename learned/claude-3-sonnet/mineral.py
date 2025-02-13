"""
Classifies: CHEBI:46662 mineral
"""
"""
Classifies: CHEBI:27008 mineral
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_mineral(smiles: str):
    """
    Determines if a molecule is a mineral based on its SMILES string.
    Minerals are typically inorganic crystalline solids formed through geological processes.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a mineral, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for organics (molecules with carbon atoms)
    if sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6) > 0:
        return False, "Contains carbon atoms, minerals are typically inorganic"

    # Check for ionic bonds (metal-nonmetal bonds)
    has_ionic_bonds = False
    for bond in mol.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        if (atom1.GetAtomicNum() in range(19, 57) and atom2.GetAtomicNum() in range(57, 87)) or \
           (atom2.GetAtomicNum() in range(19, 57) and atom1.GetAtomicNum() in range(57, 87)):
            has_ionic_bonds = True
            break
    if not has_ionic_bonds:
        return False, "No ionic bonds found, minerals typically contain metal-nonmetal bonds"

    # Check for common mineral elements
    common_mineral_elements = [8, 9, 11, 12, 13, 14, 15, 16, 17, 19, 20, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 33, 34, 35, 37, 38, 39, 40, 42, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 55, 56, 57, 58, 60, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 90, 92]
    if not any(atom.GetAtomicNum() in common_mineral_elements for atom in mol.GetAtoms()):
        return False, "Does not contain common mineral elements"

    return True, "Inorganic compound with ionic bonds and common mineral elements"