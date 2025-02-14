"""
Classifies: CHEBI:38958 indole alkaloid
"""
from rdkit import Chem

def is_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is an indole alkaloid based on its SMILES string.
    An indole alkaloid contains an indole skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an indole alkaloid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # More specific SMARTS pattern for indole:
    #  - c1: Aromatic carbon
    #  - c2: Aromatic carbon
    #  - [nH+] or [nH]: Nitrogen (protonated or not)
    #  - c3-c6: Remaining aromatic carbons
    #  - The ~ represents an aromatic bond and is more explicit than specifying double and single bonds
    indole_pattern = Chem.MolFromSmarts("c1~c~[nH+]~c2~c~c~c~c1~2") # Added protonated N
    indole_pattern2 = Chem.MolFromSmarts("c1~c~[nH]~c2~c~c~c~c1~2") # added neutral N

    if mol.HasSubstructMatch(indole_pattern) or mol.HasSubstructMatch(indole_pattern2):
        return True, "Contains an indole substructure"
    else:
        return False, "Does not contain an indole substructure"