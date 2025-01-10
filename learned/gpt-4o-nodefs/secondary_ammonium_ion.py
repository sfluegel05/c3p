"""
Classifies: CHEBI:137419 secondary ammonium ion
"""
from rdkit import Chem

def is_secondary_ammonium_ion(smiles: str):
    """
    Determines if a molecule is a secondary ammonium ion based on its SMILES string.
    A secondary ammonium ion has a positively charged nitrogen atom bonded to two organic substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary ammonium ion, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for secondary ammonium ion
    # A secondary ammonium ion is an [NH2+] (positively charged nitrogen with hydrogens)
    secondary_ammonium_pattern = Chem.MolFromSmarts("[NH2+]([C])[C]")

    # Check for matches
    if mol.HasSubstructMatch(secondary_ammonium_pattern):
        # Validate that the nitrogen atom has exactly two carbon substitutions
        matches = mol.GetSubstructMatches(secondary_ammonium_pattern)
        for match in matches:
            n_atom = mol.GetAtomWithIdx(match[0])
            bonded_carbons = [neighbor.GetAtomicNum() for neighbor in n_atom.GetNeighbors() if neighbor.GetAtomicNum() == 6]
            if len(bonded_carbons) == 2:
                return True, "Contains secondary ammonium ion with [NH2+] and two organic (carbon) substituents"

    return False, "No secondary ammonium ion structure found"