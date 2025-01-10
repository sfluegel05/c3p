"""
Classifies: CHEBI:33916 aldopentose
"""
from rdkit import Chem

def is_aldopentose(smiles: str):
    """
    Determines if a molecule is an aldopentose based on its SMILES string.

    An aldopentose is a pentose (5-carbon sugar) with a potential aldehyde
    group at one end. Commonly exists in linear or cyclic furanose/pyranose forms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an aldopentose, False otherwise.
        str: Reason for classification.
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for exactly 5 carbon atoms (pentose sugar)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 5:
        return False, f"Expected 5 carbon atoms, found {c_count}"

    # SMARTS pattern for aldehyde group (linear form)
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)[CH2]")  # Linear terminal aldehyde
    if mol.HasSubstructMatch(aldehyde_pattern):
        return True, "Contains aldehyde group, indicating a linear form of aldopentose"

    # Look for cyclic furanose (5-membered) or pyranose (6-membered) rings
    furanose_pattern = Chem.MolFromSmarts("[O][C][C][C][O]")
    pyranose_pattern = Chem.MolFromSmarts("[O][C][C][C][C][O]")
    
    if mol.HasSubstructMatch(furanose_pattern):
        return True, "Contains cyclic furanose structure"

    if mol.HasSubstructMatch(pyranose_pattern):
        return True, "Contains cyclic pyranose structure"

    # If no criteria match, it's not an aldopentose
    return False, "Does not match aldopentose structure requirements"