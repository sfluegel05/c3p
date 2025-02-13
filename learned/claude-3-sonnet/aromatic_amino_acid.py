"""
Classifies: CHEBI:33856 aromatic amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_aromatic_amino_acid(smiles: str):
    """
    Determines if a molecule is an aromatic amino acid based on its SMILES string.
    An aromatic amino acid is an amino acid with an aromatic ring in its side chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an aromatic amino acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify amino acid backbone
    amino_acid_pattern = Chem.MolFromSmarts("[N;H2,H1;!$(N-[!#6]);!$(N-[!#6]=[!#6])]C(=O)[O;H1,-]")
    if not mol.HasSubstructMatch(amino_acid_pattern):
        return False, "Molecule does not contain an amino acid backbone"

    # Check for aromatic ring attachment
    aromatic_amino_acid_pattern = Chem.MolFromSmarts("[N;H2,H1;!$(N-[!#6]);!$(N-[!#6]=[!#6])]C(=O)[O;H1,-]C1=[c;r]")
    if not mol.HasSubstructMatch(aromatic_amino_acid_pattern):
        return False, "Aromatic ring not attached to the amino acid backbone"

    # Consider stereochemistry (if needed)
    # You can use RDKit's built-in functions to determine the stereochemistry
    # and handle specific stereoisomers if required

    # Handle edge cases (e.g., multiple aromatic rings, aromatic heterocycles)
    # You can add additional checks or patterns for specific edge cases

    return True, "Molecule is an aromatic amino acid"