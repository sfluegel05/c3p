"""
Classifies: CHEBI:15889 sterol
"""
from rdkit import Chem

def is_sterol(smiles: str):
    """
    Determines if a molecule is a sterol based on its SMILES string.
    A sterol is a sterane skeleton (ABCD ring system) with a hydroxyl group at the 3 position and potential side chains.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sterol, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for sterane skeleton (ABCD ring system)
    sterane_pattern = Chem.MolFromSmarts('C1CCC2C(C1)CCC3C2CCC4(C3)CCCCC4')
    if not mol.HasSubstructMatch(sterane_pattern):
        return False, "No sterane skeleton found"
    
    # SMARTS pattern for hydroxyl group at position 3 of ring A
    hydroxyl_pattern = Chem.MolFromSmarts('[C&r3]CO')
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No 3-hydroxy group found"

    # Ensure there is additional side chain variability or functional groups
    # Validate presence of possible variable side chain
    variable_side_chain = any(atom.GetDegree() > 3 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not variable_side_chain:
        return False, "No suitable side chain found"

    return True, "Contains sterane skeleton with 3-hydroxy group and potential side chains"

# Examples given can be tested with this function to check the classification