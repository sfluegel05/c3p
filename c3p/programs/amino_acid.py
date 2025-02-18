"""
Classifies: CHEBI:33709 amino acid
"""
from rdkit import Chem

def is_amino_acid(smiles: str):
    """
    Determines if a molecule is an amino acid based on its SMILES string.
    An amino acid is characterized by having a carboxylic acid group and an amino group
    attached to the same alpha carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classifiable as an amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Carboxylic acid pattern (C(=O)O)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    carboxylic_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if not carboxylic_matches:
        return False, "No carboxylic acid group found"

    # Amino group pattern (NH2)
    amino_group_pattern = Chem.MolFromSmarts("N")
    amino_group_matches = mol.GetSubstructMatches(amino_group_pattern)
    if not amino_group_matches:
        return False, "No amino group found"
    
    # Identify if the amino group is connected to the same carbon connected to the carboxyl group
    for carbox_idx in carboxylic_matches:
        carbox_c_idx = carbox_idx[0]  # Index of the carbon in the carboxylic acid
        for amino_idx in amino_group_matches:
            amino_n_idx = amino_idx[0]  # Index of the nitrogen in the amino group
            # Check if the carbon with the carboxyl is connected directly to the nitrogen
            carbox_c_atom = mol.GetAtomWithIdx(carbox_c_idx)
            if any(nbr.GetIdx() == amino_n_idx for nbr in carbox_c_atom.GetNeighbors()):
                return True, "Contains carboxylic acid and amino group connected to same carbon (alpha carbon)"

    return False, "No amino acid backbone structure (amino and carboxyl on same C) found"