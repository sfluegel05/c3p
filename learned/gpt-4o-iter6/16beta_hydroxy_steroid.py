"""
Classifies: CHEBI:17354 16beta-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_16beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 16beta-hydroxy steroid based on its SMILES string.
    A 16beta-hydroxy steroid is defined by the presence of a hydroxyl group at position 16
    in the steroid backbone with a beta-configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 16beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for steroid backbone pattern (cyclopentanoperhydrophenanthrene structure)
    steroid_pattern = Chem.MolFromSmarts("C1CCC2C1CCC3C2CC4CC[C@]34C")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"

    # Check for the 16beta-hydroxy group
    beta_hydroxy_16_pattern = Chem.MolFromSmarts("[C@@H](O)[C@H](C)C")
    if not mol.HasSubstructMatch(beta_hydroxy_16_pattern):
        return False, "No 16beta-hydroxy group found"

    # Verify additional stereochemistry for typical steroidal structure
    stereo_matches = mol.GetSubstructMatches(beta_hydroxy_16_pattern)
    if any(mol.GetAtomWithIdx(match[1]).GetChiralTag() != Chem.CHI_TETRAHEDRAL_CCW for match in stereo_matches):
        return False, "Incorrect stereochemistry at position 16"

    return True, "Contains a steroid backbone with 16beta-hydroxy group"

# Example usage
print(is_16beta_hydroxy_steroid("C[C@]12CC[C@H]3[C@@H](CCC4=CC(=O)CC[C@]34C)[C@@H]1C[C@H](O)[C@@H]2O"))