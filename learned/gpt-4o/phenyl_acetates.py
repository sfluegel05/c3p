"""
Classifies: CHEBI:140310 phenyl acetates
"""
from rdkit import Chem

def is_phenyl_acetates(smiles: str):
    """
    Determines if a molecule is a phenyl acetate based on its SMILES string.
    A phenyl acetate is an acetate ester attached to a phenolic moiety.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is a phenyl acetate, False otherwise.
        str: Reason for classification.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the acetate ester pattern
    acetate_pattern = Chem.MolFromSmarts("C(=O)O")  # Acetate ester part

    # Define the phenol (phenyl-) moiety pattern (generic phenolic structure)
    phenol_pattern = Chem.MolFromSmarts("c1cccc(O)c1")  # Phenol pattern

    # Check if the molecule has both patterns
    if not mol.HasSubstructMatch(acetate_pattern):
        return False, "No acetate ester group found"
    
    if not mol.HasSubstructMatch(phenol_pattern):
        return False, "No phenol moiety found"

    # Ensure the structures are linked appropriately
    # Find matches for acetate and phenol patterns
    acetate_matches = mol.GetSubstructMatches(acetate_pattern)
    phenol_matches = mol.GetSubstructMatches(phenol_pattern)

    # Check connectivity: There should be a common atom where these patterns connect
    for acetate in acetate_matches:
        for phenol in phenol_matches:
            # Extract bonding atom (oxygen connection)
            acetate_oxygen = acetate[2]  # Oxygen in acetate pattern (at index 2 in C(=O)O)
            phenol_oxygen_candidates = [atom.GetIdx() for atom in mol.GetAtoms() if atom.IsInRingSize(6) and atom.GetAtomicNum() == 8]

            if acetate_oxygen in phenol_oxygen_candidates:
                return True, "Contains a phenolic moiety linked to acetate ester"

    return False, "Acetate ester and phenolic moiety are not properly connected"

# Example usage
print(is_phenyl_acetates("CC(=O)Oc1ccccc1"))  # Should return True