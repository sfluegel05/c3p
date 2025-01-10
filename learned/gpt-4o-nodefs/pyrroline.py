"""
Classifies: CHEBI:23763 pyrroline
"""
from rdkit import Chem

def is_pyrroline(smiles: str):
    """
    Determines if a molecule is a pyrroline based on its SMILES string.
    A pyrroline is characterized by a five-membered ring containing a nitrogen atom and at least one C=C double bond.

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

    # Define a flexible pyrroline SMARTS pattern: five-membered ring with a nitrogen and a double bond
    pyrroline_pattern = Chem.MolFromSmarts("[C,c;!R][C,c;!R][N,n;!R][C,c;!R]=[C,c;!R]")

    # Check if the molecule contains the generalized pyrroline substructure
    if mol.HasSubstructMatch(pyrroline_pattern):
        return True, "Contains a pyrroline structure: five-membered ring with a nitrogen and a double bond"

    return False, "Lacks pyrroline structure: specific five-membered ring criteria not satisfied"