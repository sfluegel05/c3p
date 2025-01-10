"""
Classifies: CHEBI:32876 tertiary amine
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_tertiary_amine(smiles: str):
    """
    Determines if a molecule is a tertiary amine based on its SMILES string.
    A tertiary amine is characterized by a nitrogen atom bonded to exactly 
    three carbon atoms, considering potential aromatic and aliphatic contexts.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tertiary amine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"
    
    # SMARTS pattern for tertiary amine - N bonded to three carbon atoms
    tertiary_amine_pattern = Chem.MolFromSmarts("[NX3]([C])[C][C]")
    
    # Check if the molecule matches the tertiary amine pattern
    if mol.HasSubstructMatch(tertiary_amine_pattern):
        return True, "Molecule contains a tertiary amine (N bonded to 3 carbon atoms)"
    else:
        return False, "No tertiary amine found (N not bonded to 3 carbon atoms)"

    return None, "Function execution did not reach expected returns"