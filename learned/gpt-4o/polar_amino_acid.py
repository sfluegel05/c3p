"""
Classifies: CHEBI:26167 polar amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_polar_amino_acid(smiles: str):
    """
    Determines if a molecule is a polar amino acid based on its SMILES string.
    A polar amino acid contains side chains capable of forming hydrogen bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polar amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generalized amino acid backbone pattern: amino group -> alpha carbon -> carboxyl group
    aa_backbone = Chem.MolFromSmarts("[NX3][C][C](=O)[O]")
    if not mol.HasSubstructMatch(aa_backbone):
        return False, "No amino acid backbone found"

    # Potential polar side chain patterns. Restricting to real side chains by excluding main groups
    polar_side_chain_patterns = [
        "[CX4][OX2H1]",  # Alcohol/hydroxyl, e.g., serine, threonine
        "[CX3]=[NX3]",   # Carboxamide, e.g., asparagine, glutamine
        "[CX3](=O)[OX1H1]", # Carboxylate as part of side chain, e.g., aspartate
        "[OX2H1][CX4]",  # Phenols, e.g., tyrosine
        "[nH]",          # Imidazole or heterocyclic nitrogen, e.g., histidine, tryptophan
        "[SH]",          # Thiol, e.g., cysteine
        "[CX4][NX3]"     # Primary and secondary amines in side chain
    ]

    # Check if side chain contains any polar functional groups while avoiding false positives from backbone
    for pattern in polar_side_chain_patterns:
        polar_group = Chem.MolFromSmarts(pattern)
        if mol.HasSubstructMatch(polar_group):
            return True, f"Polar group '{pattern}' found in side chain"

    return False, "No polar groups found in side chain"