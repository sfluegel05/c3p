"""
Classifies: CHEBI:61051 lipid hydroperoxide
"""
from rdkit import Chem

def is_lipid_hydroperoxide(smiles: str):
    """
    Determines if a molecule is a lipid hydroperoxide based on its SMILES string.
    A lipid hydroperoxide is a lipid molecule with one or more hydroperoxy (-OOH) groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lipid hydroperoxide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for hydroperoxy group pattern (-O-O-H)
    hydroperoxy_pattern = Chem.MolFromSmarts("OO")
    if not mol.HasSubstructMatch(hydroperoxy_pattern):
        return False, "No hydroperoxy group (-OOH) found"

    # Look for long carbon chain or unsaturated lipid-like structure
    # Commonly seen as extended chains with unsaturations
    long_chain_pattern = Chem.MolFromSmarts("CC(C)C=CCC=CCC=CC")
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "No long carbon chain or suitable lipid-like structure found"

    # Check if hydroperoxy group is attached to the carbon backbone
    for match in mol.GetSubstructMatches(hydroperoxy_pattern):
        # Ensure it connects to a carbon chain, focused on chain connectivity
        if all(atom.GetAtomicNum() == 6 for atom in mol.GetAtomWithIdx(match[0]).GetNeighbors()):
            continue
        return False, "Hydroperoxy group not attached correctly to carbon backbone"
    
    return True, "Contains hydroperoxy groups attached to a lipid backbone"