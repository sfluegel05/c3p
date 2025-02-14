"""
Classifies: CHEBI:21575 N-acetyl-amino acid
"""
"""
Classifies: N-acetyl amino acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_N_acetyl_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acetyl amino acid based on its SMILES string.
    An N-acetyl amino acid is an amino acid where the nitrogen atom is acetylated (has an acetyl group attached).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an N-acetyl amino acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS pattern for N-acetylated nitrogen
    # This matches a nitrogen atom connected to a carbonyl carbon which is connected to a methyl group (acetyl group)
    n_acetyl_pattern = Chem.MolFromSmarts("N[C](=O)C")
    if not mol.HasSubstructMatch(n_acetyl_pattern):
        return False, "No N-acetylated nitrogen found"

    # SMARTS pattern for amino acid backbone
    # This matches an alpha carbon (chiral or achiral) connected to nitrogen and a carboxyl group
    amino_acid_backbone_pattern = Chem.MolFromSmarts("[NX3][C@@H](*)C(=O)O")
    if not mol.HasSubstructMatch(amino_acid_backbone_pattern):
        return False, "No amino acid backbone found"

    # Ensure that the nitrogen with the N-acetyl group is the same nitrogen in the amino acid backbone
    n_acetyl_matches = mol.GetSubstructMatches(n_acetyl_pattern)
    amino_acid_matches = mol.GetSubstructMatches(amino_acid_backbone_pattern)

    # Check if the nitrogen atom indices overlap
    n_acetyl_nitrogens = [match[0] for match in n_acetyl_matches]
    amino_acid_nitrogens = [match[0] for match in amino_acid_matches]
    common_nitrogens = set(n_acetyl_nitrogens) & set(amino_acid_nitrogens)
    if not common_nitrogens:
        return False, "N-acetyl group is not attached to the amino acid nitrogen"

    return True, "Molecule is an N-acetyl amino acid"