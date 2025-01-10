"""
Classifies: CHEBI:33704 alpha-amino acid
"""
"""
Classifies: CHEBI:33704 alpha-amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid based on its SMILES string.
    An alpha-amino acid has an amino group (-NH2) and a carboxyl group (-COOH) attached to the same carbon atom (alpha carbon).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the alpha-amino acid pattern: [C]([NH2])([C](=O)O)
    alpha_amino_acid_pattern = Chem.MolFromSmarts("[CX4]([NH2X3])([CX3](=[OX1])[OX2H1])")
    if not mol.HasSubstructMatch(alpha_amino_acid_pattern):
        return False, "No alpha-amino acid pattern found (NH2 and COOH on the same carbon)"

    # Check for the presence of at least one amino group (-NH2)
    amino_pattern = Chem.MolFromSmarts("[NH2X3]")
    if not mol.HasSubstructMatch(amino_pattern):
        return False, "No amino group (-NH2) found"

    # Check for the presence of at least one carboxyl group (-COOH)
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxyl group (-COOH) found"

    # Ensure the amino and carboxyl groups are attached to the same carbon (alpha carbon)
    alpha_carbon_pattern = Chem.MolFromSmarts("[CX4]([NH2X3])([CX3](=[OX1])[OX2H1])")
    if not mol.HasSubstructMatch(alpha_carbon_pattern):
        return False, "Amino and carboxyl groups are not attached to the same carbon (alpha carbon)"

    return True, "Contains an amino group (-NH2) and a carboxyl group (-COOH) attached to the same carbon (alpha carbon)"