"""
Classifies: CHEBI:32876 tertiary amine
"""
from rdkit import Chem

def is_tertiary_amine(smiles: str):
    """
    Determines if a molecule is a tertiary amine based on its SMILES string.
    A tertiary amine is defined as a nitrogen atom connected to three hydrocarbyl groups, which in SMILES terms means a nitrogen bonded to three carbon atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a tertiary amine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string to obtain a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS pattern for tertiary amine: Nitrogen with 3 carbons not being part of carbonyl (C=O) or cyano groups (C#N)
    tertiary_amine_pattern = Chem.MolFromSmarts("[$([NX3](C)(C)C)]!@[$(C=O)]!@[#7]")

    if mol.HasSubstructMatch(tertiary_amine_pattern):
        return True, "Nitrogen bonded to three carbon atoms found indicating a tertiary amine"

    return False, "No nitrogen bonded to three carbon atoms found or specific configurations excluded"