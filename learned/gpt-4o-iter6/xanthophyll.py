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
    chain_pattern = Chem.MolFromSmarts("C=C" + "([#6]=[#6])" * 8)  # Looser pattern for flexibility
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "No long conjugated chain pattern typical for carotenoids"

    # Check for presence of oxygen atoms (oxygenation)
    o_pattern = Chem.MolFromSmarts("[OX2H1,OX1,C=O]")  # Detects -OH, =O, and ethers
    if not mol.HasSubstructMatch(o_pattern):
        return False, "Carotene derivative must be oxygenated"

    # Check for cyclic structures typically found in many xanthophylls
    ring_matches = [ring for ring in mol.GetRingInfo().AtomRings() if all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring)]
    if len(ring_matches) < 1:
        return False, "Carotenoid backbone typically includes cyclic rings"

    return True, "Contains features of an oxygenated carotenoid (xanthophyll)"

# Example usage:
# result, reason = is_xanthophyll("C\C(\C=C\C=C(/C)\C=C\C1=C(C)C[C@@H](O)CC1(C)C)=C/C=C/C=C(\C)/C=C/C=C(\C)/C=C/C1=C(C)C(=O)CCC1(C)C")
# print(result, reason)