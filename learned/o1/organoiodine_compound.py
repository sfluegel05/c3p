"""
Classifies: CHEBI:37142 organoiodine compound
"""
"""
Classifies: organoiodine compound
"""

from rdkit import Chem
from rdkit.Chem.rdchem import BondType

def is_organoiodine_compound(smiles: str):
    """
    Determines if a molecule is an organoiodine compound based on its SMILES string.
    An organoiodine compound contains at least one carbon-iodine (C-I) bond.

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

    # Iterate through atoms to find iodine atoms
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 53:  # Iodine atom
            formal_charge = atom.GetFormalCharge()
            degree = atom.GetDegree()  # Number of connected atoms
            # Standard iodine in organoiodine compounds should have no charge and degree of 1
            if formal_charge == 0 and degree == 1:
                neighbor = atom.GetNeighbors()[0]
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                # Check if neighbor is carbon and bond is single
                if neighbor.GetAtomicNum() == 6 and bond.GetBondType() == BondType.SINGLE:
                    return True, "Contains carbon-iodine bond with iodine in standard oxidation state (-1)"
    return False, "No suitable carbon-iodine bonds found"