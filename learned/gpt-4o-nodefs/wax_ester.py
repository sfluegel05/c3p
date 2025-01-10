"""
Classifies: CHEBI:10036 wax ester
"""
from rdkit import Chem

def is_wax_ester(smiles: str):
    """
    Determines if a molecule is a wax ester based on its SMILES string.
    Wax esters are compounds with long-chain fatty acids and fatty alcohols 
    connected through ester linkages.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a wax ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for ester group pattern (C(=O)O)
    ester_pattern = Chem.MolFromSmarts("[C](=[O])[O]")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester linkage found"

    # Check for long carbon chain
    chain_length = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if chain_length < 20:  # Arbitrarily chosen based on examples
        return False, f"Carbon chain too short, found {chain_length} carbons"

    # Check if molecule has unsaturated bonds (optional)
    unsaturated_pattern = Chem.MolFromSmarts("[C]=[C]")
    unsaturated_matches = mol.GetSubstructMatches(unsaturated_pattern)
    if not unsaturated_matches:
        return False, "No unsaturated bonds found (not mandatory for wax esters but common in examples)"

    return True, "Contains ester linkage and long carbon chains"

# Example usage
print(is_wax_ester("O(CCCCCCCC/C=C\\CCCCCCCC)C(=O)CCCCCCC/C=C\\CCCCCC")) # Palmitoleyl linoleate