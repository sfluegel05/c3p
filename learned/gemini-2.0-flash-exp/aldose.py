"""
Classifies: CHEBI:15693 aldose
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_aldose(smiles: str):
    """
    Determines if a molecule is an aldose based on its SMILES string.
    Aldoses are polyhydroxy aldehydes (H[CH(OH)]nC(=O)H, n >= 2) and their intramolecular hemiacetals.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aldose, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbon atoms
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if len(carbon_atoms) < 3:
        return False, "Aldoses must contain at least 3 carbon atoms"

    # Check for aldehyde group
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)")
    has_aldehyde = mol.HasSubstructMatch(aldehyde_pattern)

    # Check for intramolecular hemiacetal
    hemiacetal_pattern = Chem.MolFromSmarts("[OX2]1[CX4H]([OX2H])[CX4H]~[CX4H]~[CX4H]1")
    has_hemiacetal = mol.HasSubstructMatch(hemiacetal_pattern)

    # Verify polyhydroxy nature and chain length
    hydroxy_pattern = Chem.MolFromSmarts("[CX4H][OX2H]")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    num_hydroxy = len(hydroxy_matches)

    if has_aldehyde:
        if num_hydroxy >= 2:
            return True, "Open-chain aldose structure detected"
        else:
            return False, "Open-chain aldehyde detected, but not enough OH groups"


    if has_hemiacetal:
       if num_hydroxy >= 3: # at least 3 for the cyclic version, one for the hemiacetal.
          return True, "Cyclic hemiacetal aldose structure detected"
       else:
           return False, "Cyclic hemiacetal detected, but not enough OH groups"


    return False, "Does not match aldose criteria"