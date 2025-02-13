"""
Classifies: CHEBI:15705 L-alpha-amino acid
"""
"""
Classifies: L-alpha-amino acid
"""
from rdkit import Chem

def is_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an L-alpha-amino acid based on its SMILES string.
    An L-alpha-amino acid is any alpha-amino acid having L-configuration at the alpha carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an L-alpha-amino acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Assign stereochemistry to the molecule
    Chem.AssignAtomChiralTagsFromStructure(mol)
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)

    # Define SMARTS pattern for L-alpha-amino acid backbone
    # [N;H2,H1]: Amino group with one or two hydrogens
    # [C@@H]: Chiral carbon with L-configuration
    # C(=O)O: Carboxylic acid group
    pattern_L = Chem.MolFromSmarts('[N;H2,H1]-[C@@H]-C(=O)O')

    # Check if the molecule matches the L-alpha-amino acid pattern
    if mol.HasSubstructMatch(pattern_L):
        return True, "Molecule matches L-alpha-amino acid backbone"
    else:
        return False, "Molecule does not match L-alpha-amino acid backbone"