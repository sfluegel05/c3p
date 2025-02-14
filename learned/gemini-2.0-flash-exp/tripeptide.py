"""
Classifies: CHEBI:47923 tripeptide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_tripeptide(smiles: str):
    """
    Determines if a molecule is a tripeptide based on its SMILES string.
    A tripeptide is a peptide consisting of three amino acid residues connected by two peptide bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tripeptide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for two peptide bonds - modified to look for peptide bonds irrespective of their termini
    peptide_bond_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[NX2]")
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    if len(peptide_bond_matches) != 2:
         return False, f"Found {len(peptide_bond_matches)} peptide bonds, expected 2"
    
    # check for free carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if len(carboxylic_acid_matches) < 1:
        return False, f"Found {len(carboxylic_acid_matches)} free carboxyl groups, expected at least 1"
    
    return True, "Contains three amino acid residues connected by two peptide bonds."