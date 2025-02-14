"""
Classifies: CHEBI:35267 quaternary ammonium ion
"""
"""
Classifies: CHEBI:35279 quaternary ammonium ion
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_quaternary_ammonium_ion(smiles: str):
    """
    Determines if a molecule is a quaternary ammonium ion based on its SMILES string.
    A quaternary ammonium ion is a nitrogen atom with a +1 charge, bonded to four substituents.

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

    # Flag to check if quaternary ammonium ion is found
    qa_found = False

    # Iterate over all atoms in the molecule
    for atom in mol.GetAtoms():
        # Check if atom is nitrogen
        if atom.GetAtomicNum() == 7:
            # Check if nitrogen has a formal charge of +1
            if atom.GetFormalCharge() == 1:
                # Check if nitrogen is bonded to four substituents
                if atom.GetDegree() == 4:
                    # Check if all bonds are single bonds
                    bond_types = [bond.GetBondType() for bond in atom.GetBonds()]
                    if all(bond_type == rdchem.BondType.SINGLE for bond_type in bond_types):
                        qa_found = True
                        break

    if qa_found:
        return True, "Contains a nitrogen atom with +1 charge bonded to four substituents (quaternary ammonium ion)"
    else:
        return False, "No quaternary ammonium ion found in the molecule"