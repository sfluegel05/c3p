"""
Classifies: CHEBI:18379 nitrile
"""
"""
Classifies: CHEBI:46419 nitrile
A compound having the structure RC#N; thus a C-substituted derivative of hydrocyanic acid, HC#N. 
In systematic nomenclature, the suffix nitrile denotes the triply bound #N atom, not the carbon atom attached to it.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_nitrile(smiles: str):
    """
    Determines if a molecule is a nitrile based on its SMILES string.
    A nitrile contains a cyano group (-C#N) attached to a carbon atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nitrile, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for cyano group (-C#N) attached to a carbon atom
    cyano_pattern = Chem.MolFromSmarts("[C][C#N]")
    if not mol.HasSubstructMatch(cyano_pattern):
        return False, "No cyano group (-C#N) attached to a carbon atom found"
    
    # Check for common non-nitrile functional groups
    amide_pattern = Chem.MolFromSmarts("[C](=O)[N]")
    if mol.HasSubstructMatch(amide_pattern):
        return False, "Contains an amide group, which is not a nitrile"
    
    ester_pattern = Chem.MolFromSmarts("[C](=O)[O]")
    if mol.HasSubstructMatch(ester_pattern):
        return False, "Contains an ester group, which is not a nitrile"
    
    # If cyano group attached to carbon and no common non-nitrile groups, classify as nitrile
    return True, "Contains a cyano group (-C#N) attached to a carbon atom"