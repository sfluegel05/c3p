"""
Classifies: CHEBI:28863 flavanones
"""
"""
Classifies: Flavanones

Determines if a molecule is a flavanone based on its SMILES string.
A flavanone has a 3,4-dihydro-2-aryl-2H-1-benzopyran-4-one skeleton.
"""

from rdkit import Chem

def is_flavanones(smiles: str):
    """
    Determines if a molecule is a flavanone based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavanone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the chroman-4-one core SMARTS pattern
    chromanone_smarts = '[O]=C1CCCOc2ccccc12'  # Chroman-4-one core fused with an aromatic ring
    chromanone_pattern = Chem.MolFromSmarts(chromanone_smarts)
    if chromanone_pattern is None:
        return None, "Error in defining chroman-4-one SMARTS pattern"

    if not mol.HasSubstructMatch(chromanone_pattern):
        return False, "Does not contain chroman-4-one core"

    # Now check if the carbon at position 2 is connected to an aryl group
    matches = mol.GetSubstructMatches(chromanone_pattern)
    for match in matches:
        # Atom indices in chromanone_smarts:
        # 0: O of ketone
        # 1: C (ketone carbon)
        # 2: C at position 2
        # 3: C
        # 4: C
        # 5: O in ring
        # 6-11: Aromatic ring fused to heterocycle

        atom_idx_C2 = match[2]  # Carbon at position 2
        atom_C2 = mol.GetAtomWithIdx(atom_idx_C2)
        for bond in atom_C2.GetBonds():
            nbr_atom = bond.GetOtherAtom(atom_C2)
            nbr_idx = nbr_atom.GetIdx()
            if nbr_idx in match:
                continue  # Skip atoms in the chromanone core
            if nbr_atom.GetIsAromatic() and nbr_atom.IsInRing():
                return True, "Contains flavanone core with aryl group at position 2"

    return False, "Does not have aryl group attached at position 2 of chromanone core"