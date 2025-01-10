"""
Classifies: CHEBI:10036 wax ester
"""
from rdkit import Chem

def is_wax_ester(smiles: str):
    """
    Determines if a molecule is a wax ester based on its SMILES string.
    Wax esters consist of long-chain fatty acids and fatty alcohols 
    connected through ester linkages, typically featuring one significant ester linkage.

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

    # Check for the presence of ester linkage (C(=O)O)
    ester_pattern = Chem.MolFromSmarts("[C](=[O])[O]")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester linkage found"

    # Estimate the number of ester linkages
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester linkages, exactly one main ester linkage expected"

    # Check for two long carbon chains on both sides of the ester group
    def count_carbon_chains(side_mol):
        return sum(1 for atom in side_mol.GetAtoms() if atom.GetAtomicNum() == 6)

    carbon_count = count_carbon_chains(mol)
    if carbon_count < 20:
        return False, "Overall carbon chain count too short"

    # Check if molecule has any unsaturated bonds and mark as a feature
    unsaturated_pattern = Chem.MolFromSmarts("[C]=[C]")
    unsaturated_matches = mol.GetSubstructMatches(unsaturated_pattern)
    unsaturated_info = "has unsaturated bonds" if unsaturated_matches else "all bonds saturated (common, not mandatory)"

    return True, f"Contains ester linkage with balanced long carbon chains; {unsaturated_info}"

# Example Usage
print(is_wax_ester("O(CCCCCCCC/C=C\\CCCCCCCC)C(=O)CCCCCCC/C=C\\CCCCCC")) # Example SMILES from the task