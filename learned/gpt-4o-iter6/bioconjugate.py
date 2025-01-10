"""
Classifies: CHEBI:64985 bioconjugate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_bioconjugate(smiles: str):
    """
    Determines if a molecule is a bioconjugate based on its SMILES string.
    A bioconjugate consists of at least 2 biological molecules covalently linked together.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bioconjugate, False otherwise
        str: Reason for classification
    """

    # Parse SMILES string to a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS patterns for biological motifs or moieties
    peptide_pattern = Chem.MolFromSmarts("N[C@@H](C)C(=O)O")   # General peptide bond
    nucleoside_pattern = Chem.MolFromSmarts("n1cnc2c1ncnc2N")  # Purine nucleobase
    cofactor_pattern = Chem.MolFromSmarts("P(=O)(O)OCCN")      # CoA-like segment

    # Check for presence of biological motifs
    peptide_matches = mol.GetSubstructMatches(peptide_pattern)
    nucleoside_matches = mol.GetSubstructMatches(nucleoside_pattern)
    cofactor_matches = mol.GetSubstructMatches(cofactor_pattern)

    # Assume we need at least two matches of any type to suggest bioconjugation
    total_biomotif_matches = len(peptide_matches) + len(nucleoside_matches) + len(cofactor_matches)
    if total_biomotif_matches < 2:
        return False, "Less than two distinct biological motifs found"

    # Check for common link groups like esters and amides
    ester_amide_pattern = Chem.MolFromSmarts("C(=O)O")  # Simplified ester/amide linkage
    ester_amide_matches = mol.HasSubstructMatch(ester_amide_pattern)
    
    if not ester_amide_matches:
        return False, "No common ester/amide linkages found"

    # Consider molecular size as an indirect measure of complexity/fusion of entities
    mol_weight = AllChem.CalcExactMolWt(mol)
    if mol_weight < 500:
        return False, "Molecular weight too low for likely bioconjugate"

    return True, "Contains multiple distinct biological motifs and linkage features, indicative of a bioconjugate"