"""
Classifies: CHEBI:32876 tertiary amine
"""
from rdkit import Chem

def is_tertiary_amine(smiles: str):
    """
    Determines if a molecule is a tertiary amine based on its SMILES string.
    A tertiary amine has a nitrogen atom bonded only through single bonds to exactly 
    three carbon atoms.

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
    
    # SMARTS pattern for tertiary amine: Nitrogen bonded to three carbon atoms, all single bonds
    tertiary_amine_pattern = Chem.MolFromSmarts("[NX3;R0][C;R0]([R0])[C;R0][C;R0]")
    
    if mol.HasSubstructMatch(tertiary_amine_pattern):
        return True, "Molecule contains a tertiary amine (N bonded to 3 carbon atoms)"
    
    return False, "No tertiary amine structure found"