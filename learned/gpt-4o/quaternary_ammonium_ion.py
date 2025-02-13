"""
Classifies: CHEBI:35267 quaternary ammonium ion
"""
from rdkit import Chem

def is_quaternary_ammonium_ion(smiles: str):
    """
    Determines if a molecule is a quaternary ammonium ion based on its SMILES string.
    A quaternary ammonium ion is defined by a nitrogen atom with a formal positive charge,
    bonded to four univalent groups, typically carbon-based.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a quaternary ammonium ion, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for quaternary ammonium feature: [N+] with four single bonds
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:  # Check if it's Nitrogen
            if atom.GetFormalCharge() == 1:  # Check if it has a positive charge
                bonded_atoms = [bond.GetOtherAtom(atom) for bond in atom.GetBonds()]
                if len(bonded_atoms) == 4:  # Ensure there are 4 bonds
                    if all(bonded_atom.GetAtomicNum() == 6 for bonded_atom in bonded_atoms):  # All to carbon
                        # Quaternary texture further validated via simple check
                        return True, "Contains a quaternary ammonium nitrogen: [N+](C)(C)(C)(C)"
                # Include any structural notes if needed in further checks
    return False, "Does not contain a quaternary ammonium nitrogen"