"""
Classifies: CHEBI:48030 tetrapeptide
"""
from rdkit import Chem

def is_tetrapeptide(smiles: str):
    """
    Determines if a molecule is a tetrapeptide based on its SMILES string.
    A tetrapeptide contains four amino-acid residues connected by peptide linkages.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrapeptide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for peptide bond pattern (C(=O)-N)
    peptide_bond_pattern = Chem.MolFromSmarts("C(=O)N")
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    
    if len(peptide_bond_matches) != 3:
        return False, f"Found {len(peptide_bond_matches)} peptide bonds, need exactly 3 for tetrapeptide"

    # Each peptide bond connects two amino-acids, total should be four amino-acids
    num_amino_acids = 4

    # Check for count of amine and carboxylic acid groups which indicates amino acids
    amine_pattern = Chem.MolFromSmarts("[NX3][CX4][CX3](=O)")
    amine_matches = mol.GetSubstructMatches(amine_pattern)

    carboxy_pattern = Chem.MolFromSmarts("C(=O)[OX2H1]")
    carboxy_matches = mol.GetSubstructMatches(carboxy_pattern)
    
    # Validate number of amino acids based on amine and carboxylate groups
    if len(amine_matches) != num_amino_acids or len(carboxy_matches) != num_amino_acids:
        return False, "Incorrect number of amine or carboxyl groups for tetrapeptide structure"
    
    return True, "Contains four amino-acid residues connected by three peptide bonds, valid tetrapeptide structure"