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

    # Identify the peptide bond pattern - carbonyl (C=O) connected to nitrogen (N)
    peptide_bond_pattern = Chem.MolFromSmarts("C(=O)N")
    peptide_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    
    # Check for exactly one central peptide (amide) bond
    if len(peptide_matches) != 1:
        return False, f"Found {len(peptide_matches)} peptide bonds, need exactly 1 for a dipeptide"

    # Check for amino acid residues: presence of at least one primary amine and one carboxyl group
    # These can be modified in the dipeptide head/tail, so allow varying degrees
    amine_pattern = Chem.MolFromSmarts("[NX3][CX3H2,CX4H3,CX3H1]=,=,,”)
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    amine_matches = len(mol.GetSubstructMatches(amine_pattern))
    carboxyl_matches = len(mol.GetSubstructMatches(carboxyl_pattern))

    # If appropriate functional groups exist to form two amino acids after accounting for bond formation
    if amine_matches < 1 or carboxyl_matches < 1:
        return False, f"Found {amine_matches} amine and {carboxyl_matches} carboxyl groups, insufficient for two amino acids"

    return True, "Contains two amino-acid residues connected by a single peptide linkage"

# This function is meant to classify dipeptides based on the structural requirements.