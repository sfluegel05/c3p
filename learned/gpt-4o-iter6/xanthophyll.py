"""
Classifies: CHEBI:27325 xanthophyll
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_xanthophyll(smiles: str):
    """
    Determines if a molecule is a xanthophyll based on its SMILES string.
    A xanthophyll is a carotenoid that is oxygenated.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a xanthophyll, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for long conjugated chain pattern
    # Long chains have alternating single and double bonds
    chain_pattern = Chem.MolFromSmarts("C=C" + ("C=C" * 10))  # At least 20 carbons with conjugating double bonds
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "No long conjugated chain pattern typical for carotenoids"

    # Check for presence of oxygen atoms
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 1:
        return False, "Carotene derivative must be oxygenated"

    # Check for the presence of cyclic ring structures typically found in many xanthophylls
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count < 1:
        return False, "Carotenoid backbone typically includes cyclic rings"

    return True, "Contains features of an oxygenated carotenoid (xanthophyll)"

# Example usage:
# is_xanthophyll("C\C(\C=C\C=C(/C)\C=C\C1=C(C)C[C@@H](O)CC1(C)C)=C/C=C/C=C(\C)/C=C/C=C(\C)/C=C/C1=C(C)C(=O)CCC1(C)C")