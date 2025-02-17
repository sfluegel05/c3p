"""
Classifies: CHEBI:33709 amino acid
"""
"""
Classifies chemical entities of the class amino acid.
Definition: A carboxylic acid containing one or more amino groups.
Examples include (S)-gabaculine, L-thyroxine, 1-aminocyclopropanecarboxylic acid, etc.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_amino_acid(smiles: str):
    """
    Determines if a molecule is an amino acid based on its SMILES string.
    A molecule qualifies as an amino acid if it contains at least one carboxylic acid
    group (either in its neutral or deprotonated form) and one or more amino groups
    (non-amide, sp3-hybridized nitrogen).

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as an amino acid, False otherwise.
        str: Explanation for the decision.
    """
    # Try to parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS for carboxylic acid groups.
    # Neutral carboxylic acid: C(=O)[OH] (hydroxyl group present)
    carboxyl_neutral = Chem.MolFromSmarts("[CX3](=O)[OX2H]")
    # Anionic carboxylate: C(=O)[O-] (deprotonated form)
    carboxyl_anion = Chem.MolFromSmarts("[CX3](=O)[O-]")
    
    # Check if either carboxyl group is present.
    has_carboxyl = mol.HasSubstructMatch(carboxyl_neutral) or mol.HasSubstructMatch(carboxyl_anion)
    if not has_carboxyl:
        return False, "No carboxylic acid group found"
    
    # Define SMARTS for amino groups.
    # We focus on sp3 nitrogen atoms (which typically represent free amino groups)
    # and try to avoid matching amide nitrogen atoms (which are usually sp2).
    amino_group = Chem.MolFromSmarts("[NX3;!$(NC=O)]")
    
    if not mol.HasSubstructMatch(amino_group):
        return False, "No amino group found"

    return True, "Found carboxylic acid group and at least one amino group"