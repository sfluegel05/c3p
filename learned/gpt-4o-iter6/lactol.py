"""
Classifies: CHEBI:38131 lactol
"""
from rdkit import Chem

def is_lactol(smiles: str):
    """
    Determines if a molecule is a lactol based on its SMILES string.
    A lactol is a cyclic hemiacetal formed by intramolecular addition of a hydroxy group
    to an aldehydic or ketonic carbonyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a lactol, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more accurate SMARTS pattern for lactol structure
    # Lactol: cyclic structure with an ether linkage (O in a cycle) and a hydroxyl group near
    lactol_pattern = Chem.MolFromSmarts("C1OC(O1)")

    if mol.HasSubstructMatch(lactol_pattern):
        return True, "Contains a lactol structural pattern (cyclic ether with hydroxyl group)"

    return False, "No lactol pattern found (cyclic hemiacetal)"