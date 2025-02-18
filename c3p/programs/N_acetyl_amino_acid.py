"""
Classifies: CHEBI:21575 N-acetyl-amino acid
"""
"""
Classifies: CHEBI:21547 N-acetyl-amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_N_acetyl_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acetyl-amino acid based on its SMILES string.
    An N-acetyl-amino acid is an N-acyl-amino acid that has acetyl as the acyl group.

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

    # Look for the acetyl group attached to a nitrogen (N-acetyl)
    n_acetyl_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[CX4H3]")
    n_acetyl_matches = mol.GetSubstructMatches(n_acetyl_pattern)
    if len(n_acetyl_matches) == 0:
        return False, "No acetyl group attached to nitrogen found"

    # Look for the amino acid backbone (carboxyl group -C(=O)O and amino group -NH-)
    amino_acid_pattern = Chem.MolFromSmarts("[NX3][CX4H][CX3](=[OX1])[OX2H1]")
    amino_acid_matches = mol.GetSubstructMatches(amino_acid_pattern)
    if len(amino_acid_matches) == 0:
        return False, "No amino acid backbone found"

    # Check if the acetyl group is directly attached to the nitrogen of the amino acid backbone
    n_acetyl_amino_acid_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[CX4H3].[NX3][CX4H][CX3](=[OX1])[OX2H1]")
    n_acetyl_amino_acid_matches = mol.GetSubstructMatches(n_acetyl_amino_acid_pattern)
    if len(n_acetyl_amino_acid_matches) == 0:
        return False, "Acetyl group not attached to the nitrogen of the amino acid backbone"

    return True, "Contains an acetyl group attached to the nitrogen of an amino acid backbone"