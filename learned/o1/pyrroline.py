"""
Classifies: CHEBI:23763 pyrroline
"""
"""
Classifies: CHEBI:25122 pyrroline
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_pyrroline(smiles: str):
    """
    Determines if a molecule is a pyrroline based on its SMILES string.
    A pyrroline is any organic heteromonocyclic compound with a structure based on dihydropyrrole,
    which is a five-membered, non-aromatic ring containing one nitrogen atom and one double bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a pyrroline, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for pyrroline
    # Non-aromatic five-membered ring with one nitrogen atom and one double bond
    pyrroline_smarts = '[#7;$([n&H0]),$([n&H1]);R5;$([R]=1)]@[#6;R1]@[#6;R1]=[#6;R1]@[#6;R1]@1'

    pyrroline_pattern = Chem.MolFromSmarts(pyrroline_smarts)
    if pyrroline_pattern is None:
        return False, "Invalid SMARTS pattern"

    # Search for pyrroline substructure
    if mol.HasSubstructMatch(pyrroline_pattern):
        return True, "Contains a pyrroline ring"

    # Alternatively, check for all possible pyrroline isomers
    # SMARTS patterns for 2-pyrroline and 3-pyrroline
    pyrroline_smarts_list = [
        'C1=CCCN1',  # 2-Pyrroline
        'C1CC=CN1',  # 3-Pyrroline
        'C1=CC=CN1', # Not a pyrroline (pyrrole), excluded
    ]

    for smarts in pyrroline_smarts_list:
        pattern = Chem.MolFromSmarts(smarts)
        if pattern and mol.HasSubstructMatch(pattern):
            return True, "Contains a pyrroline ring"

    return False, "Does not contain a pyrroline ring"