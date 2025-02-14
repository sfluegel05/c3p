"""
Classifies: CHEBI:33704 alpha-amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid based on its SMILES string.
    An alpha-amino acid has an amino group and a carboxylic acid group attached to the same carbon atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an alpha-amino acid, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS pattern for an alpha-amino acid
    # This pattern specifies that the amino group (N) and carboxylic acid group (C(=O)O)
    # are attached to the SAME carbon atom (CX4) and it does not match other similar substructures.
    alpha_amino_acid_pattern = Chem.MolFromSmarts("[NX3;H0,H1,H2][CX4]([CX3](=[OX1,OX2])[OX1,OX2])")

    if not mol.HasSubstructMatch(alpha_amino_acid_pattern):
        return False, "No alpha-amino acid substructure found."
    
    return True, "Molecule matches the alpha-amino acid definition."