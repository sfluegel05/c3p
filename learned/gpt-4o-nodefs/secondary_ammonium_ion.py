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
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Define SMARTS pattern for secondary ammonium ion
    # N+ represents a positively charged nitrogen having two carbonaceous substituents
    secondary_ammonium_pattern = Chem.MolFromSmarts("[N+;H2]([C])[C]")  # Positively charged nitrogen with two organic (carbon) bonds

    # Check for matches
    if mol.HasSubstructMatch(secondary_ammonium_pattern):
        # Ensure both substituents are carbon-connected (organic)
        matches = mol.GetSubstructMatches(secondary_ammonium_pattern)
        for match in matches:
            n_atom = mol.GetAtomWithIdx(match[0])
            non_h_atom_bonds = [nbr.GetAtomicNum() for nbr in n_atom.GetNeighbors() if nbr.GetAtomicNum() != 1]
            if len(non_h_atom_bonds) == 2 and all(atom_num == 6 for atom_num in non_h_atom_bonds):  # Check if both are carbon
                return True, "Contains secondary ammonium ion with [NH2+] and two organic (carbon) substituents"
    
    return False, "No secondary ammonium ion structure found"

# Example usage:
# print(is_secondary_ammonium_ion("CC[NH2+]CC"))