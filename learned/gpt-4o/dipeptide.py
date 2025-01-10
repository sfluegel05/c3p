"""
Classifies: CHEBI:46761 dipeptide
"""
from rdkit import Chem

def is_dipeptide(smiles: str):
    """
    Determines if a molecule is a dipeptide based on its SMILES string.
    A dipeptide contains two amino-acid residues connected by a single peptide linkage (amide bond).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dipeptide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify peptide bond patterns: carbonyl connected to nitrogen
    peptide_bond_pattern = Chem.MolFromSmarts("N[CX3](=O)C")
    peptide_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    
    # Check for exactly one peptide (amide) bond
    if len(peptide_matches) != 1:
        return False, f"Found {len(peptide_matches)} peptide bonds, need exactly 1 for a dipeptide"
    
    # Identify amino acid patterns: presence of amine group and carboxyl group
    amine_pattern = Chem.MolFromSmarts("[NX3H2,NX3H]")
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[OX1H]")
    amine_matches = len(mol.GetSubstructMatches(amine_pattern))
    carboxyl_matches = len(mol.GetSubstructMatches(carboxyl_pattern))
    
    # A dipeptide should have two amines and two carboxyl groups
    if amine_matches != 2 or carboxyl_matches != 2:
        return False, f"Found {amine_matches} amine and {carboxyl_matches} carboxyl groups, need 2 of each for two amino acids"
    
    return True, "Contains two amino-acid residues connected by a single peptide linkage"

# This function is meant to classify dipeptides based on the structural requirements.