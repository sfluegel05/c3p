"""
Classifies: CHEBI:21575 N-acetyl-amino acid
"""
from rdkit import Chem

def is_N_acetyl_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acetyl-amino acid based on its SMILES string.
    An N-acetyl-amino acid has an acetyl group attached to the nitrogen of an amino acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acetyl-amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for acetyl group attached to nitrogen ('CC(=O)N')
    acetyl_pattern = Chem.MolFromSmarts('CC(=O)N')
    if not mol.HasSubstructMatch(acetyl_pattern):
        return False, "No acetyl group attached to nitrogen found"
    
    # Ensure it has an amino acid structure: a backbone with a carboxylic acid group
    amino_acid_pattern = Chem.MolFromSmarts('[NX3][C;!R]C(=O)O')
    if not mol.HasSubstructMatch(amino_acid_pattern):
        return False, "No amino acid backbone identified"

    # Check if the structure seems to be an amino acid derivative
    # Example: Ensure presence of carboxylic acid and amino groups
    carboxylic_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'O')
    nitrogen_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'N')
    
    if carboxylic_count < 2:
        return False, "Insufficient oxygen for carboxylic acid group"
    if nitrogen_count < 1:
        return False, "No nitrogen atom found for amino group"

    return True, "Contains acetyl group attached to amino acid nitrogen with carboxylic acid group"


# Example use:
# result, reason = is_N_acetyl_amino_acid('CC(=O)N[C@@H](Cc1c[nH]c2ccccc12)C(O)=O')
# print(result, reason)