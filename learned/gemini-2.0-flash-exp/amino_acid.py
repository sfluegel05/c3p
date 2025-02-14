"""
Classifies: CHEBI:33709 amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_amino_acid(smiles: str):
    """
    Determines if a molecule is an amino acid based on its SMILES string.
    An amino acid is defined as a carboxylic acid containing one or more amino groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
         tuple(bool, str): True if molecule is an amino acid, False otherwise and the reason.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS for amino group (any N with 3 bonds or a quaternary N) and carboxyl (carboxylic acid and carboxylate)
    amino_pattern = Chem.MolFromSmarts("[NX3+,NX3]") 
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[OX1-,OX2H0,OX2H1]")

    if not mol.HasSubstructMatch(amino_pattern):
        return False, "No amino group found"
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"
    
    # If we reach this point, it means there is at least one amino and one carboxyl group.
    # No further checks regarding the connection are needed according to the problem definition.
    return True, "Contains at least one amino and one carboxylic acid group"