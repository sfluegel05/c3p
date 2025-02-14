"""
Classifies: CHEBI:32877 primary amine
"""
from rdkit import Chem

def is_primary_amine(smiles: str):
    """
    Determines if a molecule is a primary amine based on its SMILES string.
    A primary amine is defined as a nitrogen atom bonded to one carbon atom (or other non-carbon atoms) and two hydrogen atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a primary amine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for a primary amine (N with 1 bond to carbon and two hydrogens)
    primary_amine_pattern = Chem.MolFromSmarts("[NX3;H2][CX4]") 
    
    # Check for the substructure match
    if mol.HasSubstructMatch(primary_amine_pattern):
        return True, "Contains a primary amine"
    else:
          # Define SMARTS pattern for a primary amine (N with 1 bond to a non-carbon)
        primary_amine_pattern2 = Chem.MolFromSmarts("[NX3;H2][!#6]")
        if mol.HasSubstructMatch(primary_amine_pattern2):
             return True, "Contains a primary amine, N attached to non-Carbon"
        else:
           return False, "Does not contain a primary amine"