"""
Classifies: CHEBI:137419 secondary ammonium ion
"""
from rdkit import Chem

def is_secondary_ammonium_ion(smiles: str):
    """
    Determines if a molecule is a secondary ammonium ion based on its SMILES string.
    A secondary ammonium ion has a positively charged nitrogen atom with two organic substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary ammonium ion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Define SMARTS pattern for secondary ammonium ion (considering positive nitrogen with two carbon bonds)
    secondary_ammonium_pattern = Chem.MolFromSmarts("[N+;R0]([!#1;R0])[!#1;R0]")  # Positive nitrogen with exactly two non-hydrogen substituents

    # Check matches
    if mol.HasSubstructMatch(secondary_ammonium_pattern):
        # Double check substituents count and their nature
        matches = mol.GetSubstructMatches(secondary_ammonium_pattern)
        for match in matches:
            n_atom = mol.GetAtomWithIdx(match[0])
            num_non_hydrogen_bonds = sum(1 for neighbor in n_atom.GetNeighbors() if neighbor.GetAtomicNum() != 1)
            if num_non_hydrogen_bonds == 2:  # Confirm two organic substituents
                return True, "Contains secondary ammonium ion with [NH2+] and two organic substituents"
    
    return False, "No secondary ammonium ion structure found"

# Example usage:
# print(is_secondary_ammonium_ion("CC[NH2+]CC"))