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

    # SMARTS for amino group (primary, secondary, tertiary, protonated)
    # explicitly specify hydrogens for explicit match or include charged N
    amino_pattern = Chem.MolFromSmarts("[NX3H2,NX3H1,NX3H0+,NX3+,NX2H1,NX2H0]")
    # SMARTS for carboxylic acid group, excluding those in amides and esters
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[OX1H0,OX2H1]")
    

    if not mol.HasSubstructMatch(amino_pattern):
        return False, "No amino group found"
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"
    
    # If we reach this point, it means there is at least one amino and one carboxyl group.
    return True, "Contains at least one amino and one carboxylic acid group"