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

    # Look for hydroperoxy group pattern (-OOH)
    hydroperoxy_pattern = Chem.MolFromSmarts("OO")
    if not mol.HasSubstructMatch(hydroperoxy_pattern):
        return False, "No hydroperoxy group (-OOH) found"

    # Look for a generalized long carbon chain with possible unsaturation
    long_chain_pattern = Chem.MolFromSmarts("[#6]-[#6]-[#6]-[#6](=O)~[#8]") # This pattern covers broader lipid structures
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "No long carbon chain or suitable lipid-like structure found"

    # Ensure hydroperoxy group is connected correctly to the lipid structure
    for match in mol.GetSubstructMatches(hydroperoxy_pattern):
        neighbor_atoms = [a.GetAtomicNum() for a in mol.GetAtomWithIdx(match[0]).GetNeighbors()]
        if 6 in neighbor_atoms:  # Checking if adjacent atoms include carbon
            return True, "Contains hydroperoxy groups attached to a lipid backbone"
    
    return False, "Hydroperoxy group not attached correctly to carbon backbone"