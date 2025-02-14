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

    # Refined check for aldehyde group at the end of a carbon chain
    aldehyde_pattern = Chem.MolFromSmarts("[CH1](=O)[#6]")
    has_aldehyde = mol.HasSubstructMatch(aldehyde_pattern)

    # Refined check for intramolecular hemiacetal (5 or 6 membered ring)
    hemiacetal_pattern = Chem.MolFromSmarts("[OX2]1[CX4H]([OX2H])[CX4H]([#6])[CX4H]([#6])[CX4H]1")
    has_hemiacetal = mol.HasSubstructMatch(hemiacetal_pattern)
    hemiacetal_pattern2 = Chem.MolFromSmarts("[OX2]1[CX4H]([OX2H])[CX4H]([#6])[CX4H]([#6])[CX4H]([#6])1") # 6 membered ring
    has_hemiacetal2 = mol.HasSubstructMatch(hemiacetal_pattern2)


    # Check for polyhydroxy nature
    hydroxy_pattern = Chem.MolFromSmarts("[CX4H][OX2H]")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    num_hydroxy = len(hydroxy_matches)

    #Check that the ring contains only Carbon and Oxygen (except hemiacetalic OH)
    ring_atoms_pattern = Chem.MolFromSmarts("[OX2]1[CX4H]([OX2H])[CX4H,OX2][CX4H,OX2][CX4H,OX2]1")
    ring_atoms_pattern2 = Chem.MolFromSmarts("[OX2]1[CX4H]([OX2H])[CX4H,OX2][CX4H,OX2][CX4H,OX2][CX4H,OX2]1")

    has_correct_ring_atoms = mol.HasSubstructMatch(ring_atoms_pattern)
    has_correct_ring_atoms2 = mol.HasSubstructMatch(ring_atoms_pattern2)


    # Explicit check for the polyhydroxy chain
    chain_pattern = Chem.MolFromSmarts("[CH1](=O)[CH]([OX2H])[CH]([OX2H])")
    has_polyhydroxy_chain = mol.HasSubstructMatch(chain_pattern)


    if has_aldehyde:
        if has_polyhydroxy_chain and num_hydroxy >= 2 :
            return True, "Open-chain aldose structure detected"
        else:
             return False, "Open-chain aldehyde detected, but not enough OH groups or wrong chain pattern"

    if has_hemiacetal or has_hemiacetal2:

        if num_hydroxy >=3 and (has_correct_ring_atoms or has_correct_ring_atoms2) :
           return True, "Cyclic hemiacetal aldose structure detected"
        else:
            return False, "Cyclic hemiacetal detected, but not enough OH groups or wrong ring atoms"



    return False, "Does not match aldose criteria"