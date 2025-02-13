"""
Classifies: CHEBI:15889 sterol
"""
"""
Classifies: CHEBI:15889 sterol
"""
from rdkit import Chem

def is_sterol(smiles: str):
    """
    Determines if a molecule is a sterol based on its SMILES string.
    A sterol is any 3-hydroxy steroid whose skeleton is closely related
    to cholestan-3-ol, with possible additional carbon atoms in the side chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a sterol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Pattern for a broad steroid core (ignoring specific functional groups but including 3-hydroxy)
    steroid_core_pattern = Chem.MolFromSmarts("C1CCC2C3CCC4(C)C(O)CCC4C3CCC2C1")
    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "Does not contain a sterol core structure"

    # Ensure the presence of the 3-hydroxy group
    hydroxy_pattern = Chem.MolFromSmarts("[C@H](O)C")
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "Does not have a 3-hydroxy group"

    # Include check for potential variations and substitutions at other common positions
    # (e.g., methylation, carbon chain extensions, oxidations)
    flex_steroid_pattern = Chem.MolFromSmarts("C1C[C@@]2(CC[C@H]3C4CCC([C@H](O)C4C3)C2C)C1")
    if not mol.HasSubstructMatch(flex_steroid_pattern):
        return False, "Does not fit known variants of sterol structures"

    # Verification step to ensure it doesn't exclude falsely with steroids slightly modified
    # Cholesterol is benchmark molecular example: C27H46O
    if Chem.MolFromSmiles("C1[C@@]2([C@]3(CC[C@]4([C@]([C@@]3(CC=C2C[C@H](C1)O)[H])(CC[C@@]4([C@H](C)CCCC(C)C)[H])[H])C)[H])C").HasSubstructMatch(mol):
        return True, "Recognized within known sterol framework"

    return True, "Contains a broad 3-hydroxy steroid structure related to sterol"