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
    A nitrile contains a cyano group (-C#N).

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
    
    # Look for cyano group (-C#N)
    cyano_pattern = Chem.MolFromSmarts("[C#N]")
    if not mol.HasSubstructMatch(cyano_pattern):
        return False, "No cyano group (-C#N) found"
    
    # Count cyano groups
    cyano_matches = mol.GetSubstructMatches(cyano_pattern)
    n_cyano = len(cyano_matches)
    
    # Check for multiple cyano groups
    if n_cyano > 1:
        return True, f"Contains {n_cyano} cyano groups (-C#N)"
    
    # Check for other common functional groups
    nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
    if mol.HasSubstructMatch(nitro_pattern):
        return True, "Contains a cyano group (-C#N) and a nitro group (-NO2)"
    
    amine_pattern = Chem.MolFromSmarts("[N;H]")
    if mol.HasSubstructMatch(amine_pattern):
        return True, "Contains a cyano group (-C#N) and an amine (-NH2 or -NH-)"
    
    # If only one cyano group, return True
    return True, "Contains a single cyano group (-C#N)"