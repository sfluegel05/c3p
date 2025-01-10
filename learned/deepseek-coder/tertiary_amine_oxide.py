"""
Classifies: CHEBI:134363 tertiary amine oxide
"""
"""
Classifies: CHEBI:50995 tertiary amine oxide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_tertiary_amine_oxide(smiles: str):
    """
    Determines if a molecule is a tertiary amine oxide based on its SMILES string.
    A tertiary amine oxide is defined as an N-oxide where there are three organic groups bonded to the nitrogen atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tertiary amine oxide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the N-oxide pattern: Nitrogen with three organic groups and an oxygen with a positive charge
    n_oxide_pattern = Chem.MolFromSmarts("[N+]([O-])([#6])([#6])([#6])")
    
    # Check if the molecule matches the N-oxide pattern
    if mol.HasSubstructMatch(n_oxide_pattern):
        # Ensure that the three groups bonded to the nitrogen are organic (carbon-based)
        for match in mol.GetSubstructMatches(n_oxide_pattern):
            nitrogen_atom = mol.GetAtomWithIdx(match[0])
            # Check if all three bonded atoms are carbon
            if all(mol.GetAtomWithIdx(neighbor.GetIdx()).GetAtomicNum() == 6 for neighbor in nitrogen_atom.GetNeighbors() if neighbor.GetIdx() != match[1]):
                return True, "Contains a nitrogen atom bonded to three organic groups and an oxygen atom (N-oxide)"
        
        return False, "Nitrogen is not bonded to three organic groups"
    else:
        return False, "Does not contain a nitrogen atom bonded to three organic groups and an oxygen atom (N-oxide)"