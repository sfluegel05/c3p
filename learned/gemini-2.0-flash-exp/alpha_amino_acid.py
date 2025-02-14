"""
Classifies: CHEBI:33704 alpha-amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid based on its SMILES string.
    An alpha-amino acid has an amino group and a carboxylic acid group attached to the same carbon atom (alpha carbon).
    Handles zwitterions, substitutions and salts.

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

    # Remove any counterions
    if "." in smiles:
        parts = smiles.split(".")
        if len(parts) > 1:
            # Try using the first part for core structure; this may not always work.
            mol = Chem.MolFromSmiles(parts[0])
            if mol is None:
              return False, "Invalid SMILES string after salt removal"
            
    # More specific SMARTS pattern to target the alpha carbon.
    #   [NX3,NX4+] - amino group
    #   [CX4] - alpha carbon
    #   [CX3](=[OX1,OX2])[OX1,OX2] - carboxylic acid group
    alpha_amino_acid_pattern = Chem.MolFromSmarts("[NX3,NX4+][CX4]([CX3](=[OX1,OX2])[OX1,OX2])")

    if not mol.HasSubstructMatch(alpha_amino_acid_pattern):
        return False, "No alpha-amino acid substructure found"

    return True, "Molecule matches the alpha-amino acid definition."