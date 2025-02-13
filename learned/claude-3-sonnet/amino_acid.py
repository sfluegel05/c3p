"""
Classifies: CHEBI:33709 amino acid
"""
"""
Classifies: CHEBI:33597 amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_amino_acid(smiles: str):
    """
    Determines if a molecule is an amino acid based on its SMILES string.
    An amino acid is a carboxylic acid containing one or more amino groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for carboxylic acid group (-C(=O)O)
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"
    
    # Look for amino group (-N)
    amino_pattern = Chem.MolFromSmarts("N")
    if not mol.HasSubstructMatch(amino_pattern):
        return False, "No amino group found"
    
    # Check for at least one carbon-nitrogen bond (C-N) to ensure amino group is attached
    cn_bond_pattern = Chem.MolFromSmarts("CN")
    if not AllChem.MolToSmarts(mol).count("CN") > 0:
        return False, "Amino group not attached to carbon"
    
    # Check for additional functional groups if present
    # e.g., hydroxyl (-OH), thiol (-SH), phosphate (-OP(O)(O)=O), etc.
    
    return True, "Contains a carboxylic acid group and an amino group"