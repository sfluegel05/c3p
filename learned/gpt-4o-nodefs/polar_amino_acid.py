"""
Classifies: CHEBI:26167 polar amino acid
"""
from rdkit import Chem

def is_polar_amino_acid(smiles: str):
    """
    Determines if a molecule is a polar amino acid based on its SMILES string.
    Polar amino acids have a carboxyl and amine functional group, with side chains that can form hydrogen bonds.

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

    # Redefining the core structure pattern to account for isotopic variations
    # General pattern for alpha-amino acids: N-[C]-[CX3](=O)[O]
    core_pattern = Chem.MolFromSmarts("[NX3;H2,H1;[!$(*=*)&!$(*#*)]][C;R0][CX3](=O)[OX1]")

    if not mol.HasSubstructMatch(core_pattern):
        return False, "No alpha-amino acid core structure found"

    # Improved check for polar side chains characteristic of polar amino acids
    # Focusing on specific patterns common in known polar amino acids
    polar_side_chain_patterns = [
        Chem.MolFromSmarts("[CX3](=O)[NX3H2]"),  # Amide group (Asparagine, Glutamine)
        Chem.MolFromSmarts("[CX3](=O)[[NX3]R0]"),  # General amide with no rings
        Chem.MolFromSmarts("[OX2H]"),  # Hydroxyl group (Serine, Threonine, Tyrosine)
        Chem.MolFromSmarts("[SX2H]"),  # Thiol group (Cysteine)
        Chem.MolFromSmarts("[NX3][OX1]")  # Imidazole (Histidine)
    ]

    for pattern in polar_side_chain_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains polar side chain capable of hydrogen bonding"

    return False, "No polar side chain found capable of hydrogen bonding"