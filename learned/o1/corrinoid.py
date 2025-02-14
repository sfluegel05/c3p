"""
Classifies: CHEBI:33913 corrinoid
"""
"""
Classifies: CHEBI:33909 corrinoid
"""

from rdkit import Chem

def is_corrinoid(smiles: str):
    """
    Determines if a molecule is a corrinoid based on its SMILES string.
    A corrinoid is a derivative of the corrin nucleus, which contains four reduced or partly reduced
    pyrrole rings joined in a macrocycle by three =C- groups and one direct carbon-carbon bond
    linking alpha positions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a corrinoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the corrin nucleus more accurately as a SMARTS pattern
    # The corrin nucleus consists of a macrocycle with 4 pyrrole rings connected via three methine bridges and one direct bond
    # This is a complex structure; the following is an approximation for matching purposes
    corrin_smarts = '''
    [n]1ccc2c1[n]ccc3c2c[n]cc4c3[n]ccc14
    '''
    corrin_pattern = Chem.MolFromSmarts(corrin_smarts)
    if corrin_pattern is None:
        return False, "Invalid corrin nucleus SMARTS pattern"

    # Check if the molecule contains the corrin nucleus
    if not mol.HasSubstructMatch(corrin_pattern):
        return False, "Molecule does not contain the corrin nucleus"

    # Check for the presence of a cobalt ion (atomic number 27)
    cobalt_present = any(atom.GetAtomicNum() == 27 for atom in mol.GetAtoms())
    if not cobalt_present:
        return False, "Molecule does not contain cobalt, which is characteristic of corrinoids"

    # Optionally, check that cobalt is coordinated within the corrin ring
    # This would require mapping the substructure match and verifying cobalt's position
    # For simplicity, we'll assume that the presence of cobalt and the corrin ring is sufficient

    return True, "Molecule contains the corrin nucleus and cobalt ion characteristic of corrinoids"