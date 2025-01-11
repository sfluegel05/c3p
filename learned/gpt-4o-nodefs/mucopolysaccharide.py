"""
Classifies: CHEBI:37395 mucopolysaccharide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_mucopolysaccharide(smiles: str):
    """
    Attempts to classify a molecule as a 'mucopolysaccharide' based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: Returns True if identified as mucopolysaccharide-like structure, False otherwise.
        str: Provides a reason for the classification.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Example feature check: Check for a specific substructure that can be associated 
    # with 'mucopolysaccharide-like' characteristics
    putative_features = [
        "[CX3](=O)[NX3][CX3](=O)N",  # More generalized peptide-like connections
        "O[C@H1]C=O",                # Aldehydic or ketogenic secondary features
        "O=C",                       # Carboxylate-like features
        "[NX2]C=O",                  # Amide and connected carbonyl groups
        "C[NX3]=O"                   # Amine oxides or similar
    ]
    
    # Dummy check for aforementioned features
    for feature_smarts in putative_features:
        feature_pattern = Chem.MolFromSmarts(feature_smarts)
        if feature_pattern and mol.HasSubstructMatch(feature_pattern):
            return True, f"Contains recurring substructure: {feature_smarts}"
    
    return False, "Doesn't fit identified mucopolysaccharide or known related features"