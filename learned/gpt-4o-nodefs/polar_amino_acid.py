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
    
    # Parse the SMILES string into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the core structure for an alpha-amino acid
    core_pattern = Chem.MolFromSmarts("N[C@@H](C(=O)O)")
    if not core_pattern:
        return False, "Core pattern could not be created"
    
    if not mol.HasSubstructMatch(core_pattern):
        return False, "No alpha-amino acid core structure found"

    # Define patterns for polar side chains
    polar_side_chain_patterns = [
        Chem.MolFromSmarts("C(=O)N"),  # Amide group (Asparagine, Glutamine)
        Chem.MolFromSmarts("[OH]"),    # Hydroxyl group (Serine, Threonine, Tyrosine)
        Chem.MolFromSmarts("[SH]"),    # Thiol group (Cysteine)
        Chem.MolFromSmarts("C=C([NH])[NH]")  # Imidazole ring (Histidine)
    ]

    # Check each polar side chain pattern
    for pattern in polar_side_chain_patterns:
        if pattern and mol.HasSubstructMatch(pattern):
            return True, "Contains polar side chain capable of hydrogen bonding"

    return False, "No polar side chain found capable of hydrogen bonding"