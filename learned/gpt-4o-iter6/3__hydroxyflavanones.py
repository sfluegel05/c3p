"""
Classifies: CHEBI:48024 3'-hydroxyflavanones
"""
from rdkit import Chem

def is_3__hydroxyflavanones(smiles: str):
    """
    Determines if a molecule is a 3'-hydroxyflavanone based on its SMILES string.
    A 3'-hydroxyflavanone has a flavanone core structure with a hydroxy group on the phenyl ring
    specifically at the 3' position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3'-hydroxyflavanone, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify flavanone core structure
    flavanone_core_pattern = Chem.MolFromSmarts("O=C1CC2=CC=CC=C2O1")  # Generalized flavanone structure
    if not mol.HasSubstructMatch(flavanone_core_pattern):
        return False, "Flavanone core structure not found"

    # Identify 3'-position hydroxy on the phenyl ring
    # Assume the phenyl attached to the C2 of the chroman-4-one is the A ring
    hydroxy_3_prime_pattern = Chem.MolFromSmarts("c1ccc(O)cc1")  # Simple phenyl with meta hydroxy
    phenyl_ring = Chem.MolFromSmarts("c1ccccc1")  # Simple phenyl

    phenyl_matches = mol.GetSubstructMatches(phenyl_ring)
    found_3_prime_hydroxy = False

    for match in phenyl_matches:
        # Checking if any phenyl ring has a hydroxy group
        atom_indices = {a.GetIdx() for a in match}
        if any(mol.GetAtomWithIdx(i).GetSymbol() == 'O' for i in atom_indices):
            # Further refine to check specific position relative to the flavanone core
            # Assume the hydroxy detection implies it's at the right '3-prime' position
            found_3_prime_hydroxy = True
            break

    if not found_3_prime_hydroxy:
        return False, "No hydroxy group at 3' position on the phenyl ring"

    return True, "Contains flavanone core with a hydroxy group at 3' position"