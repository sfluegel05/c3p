"""
Classifies: CHEBI:48927 N-acyl-L-alpha-amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_N_acyl_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acyl-L-alpha-amino acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acyl-L-alpha-amino acid, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for L-alpha-amino acid core pattern
    # Look for a CH that has a carboxyl group attached (-C(=O)O) and a nitrogen attached to a carbonyl
    amino_acid_core_pattern = Chem.MolFromSmarts("[NX2][CX3](=[OX1])~[CX4][CX3](=O)[OX1]")
    if not mol.HasSubstructMatch(amino_acid_core_pattern):
        return False, "Molecule does not contain L-alpha-amino acid core."

    # Look for chiral carbon (alpha carbon)
    chiral_center_pattern = Chem.MolFromSmarts("[C@H](C(=O)O)")
    if not mol.HasSubstructMatch(chiral_center_pattern):
         return False, "Molecule does not have L configuration of chiral carbon."

    return True, "Molecule is an N-acyl-L-alpha-amino acid"