"""
Classifies: CHEBI:29256 thiol
"""
from rdkit import Chem

def is_thiol(smiles: str):
    """
    Determines if a molecule is a thiol based on its SMILES string.
    A thiol is an organosulfur compound with a thiol group (-SH) attached to a carbon atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a thiol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Enhanced thiol pattern to capture thiol groups precisely: S attached to a carbon and a hydrogen
    # This prevents false positives like thioethers (sulfur bonded to two carbons) or other non-thiol compounds.
    thiol_pattern = Chem.MolFromSmarts("[#6][SH]")  # Carbon singly bonded to sulfur with hydrogen on sulfur

    # Check if a thiol group is present
    if not mol.HasSubstructMatch(thiol_pattern):
        return False, "No thiol group found"

    # Additional check for thiolate (deprotonated thiols)
    thiolate_pattern = Chem.MolFromSmarts("[#6][S-]")  # Carbon singly bonded to sulfur with negative charge

    if mol.HasSubstructMatch(thiolate_pattern):
        return True, "Contains a thiol group (or thiolate anion) attached to a carbon atom"

    return True, "Contains a thiol group (-SH) attached to a carbon atom"