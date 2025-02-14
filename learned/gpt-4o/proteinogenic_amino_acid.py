"""
Classifies: CHEBI:83813 proteinogenic amino acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_proteinogenic_amino_acid(smiles: str):
    """
    Determines if a molecule is a proteinogenic amino acid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is a proteinogenic amino acid, False otherwise.
        str: Reason for classification.
    """
    
    # Parse SMILES to create RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carboxyl group (-C(=O)O)
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[O-]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxyl group found"

    # Check for amino group (-N- or -NH2)
    amino_pattern = Chem.MolFromSmarts("[NH2,NH]")
    if not mol.HasSubstructMatch(amino_pattern):
        return False, "No amino group found"
    
    # Check for an alpha carbon connecting both groups
    alpha_carbon_pattern = Chem.MolFromSmarts("[C@H](N)(C(=O)[O-])")  # @ for chirality
    if not mol.HasSubstructMatch(alpha_carbon_pattern):
        return False, "No alpha carbon configuration found (required stereochemistry)"

    # Allow exception for glycine which is not chiral
    glycine_pattern = Chem.MolFromSmarts("C(C(=O)O)N")
    if mol.HasSubstructMatch(glycine_pattern):
        return True, "Glycine detected, which is an achiral proteinogenic amino acid"
    
    return True, "The molecule fits the profile of a proteinogenic amino acid"