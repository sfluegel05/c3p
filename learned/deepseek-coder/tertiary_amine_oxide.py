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
        # Additional check to ensure the nitrogen is not part of a more complex structure
        # such as a nitro group or other nitrogen-oxygen functional groups
        nitro_pattern = Chem.MolFromSmarts("[N+](=O)([O-])")
        if mol.HasSubstructMatch(nitro_pattern):
            return False, "Contains a nitro group, not a tertiary amine oxide"
        
        # Check if the nitrogen is part of a ring system
        for match in mol.GetSubstructMatches(n_oxide_pattern):
            nitrogen_atom = mol.GetAtomWithIdx(match[0])
            if nitrogen_atom.IsInRing():
                return False, "Nitrogen is part of a ring system, not a typical tertiary amine oxide"
        
        return True, "Contains a nitrogen atom bonded to three organic groups and an oxygen atom (N-oxide)"
    else:
        return False, "Does not contain a nitrogen atom bonded to three organic groups and an oxygen atom (N-oxide)"