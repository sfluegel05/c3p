"""
Classifies: CHEBI:25903 peptide antibiotic
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_peptide_antibiotic(smiles: str):
    """
    Determines if a molecule is a peptide antibiotic based on its SMILES string.
    Peptide antibiotics are characterized by multiple amide bonds, amino acid structures and a larger MW.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a peptide antibiotic, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for multiple amide bonds
    amide_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[NX2,NX1]")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if len(amide_matches) < 5:  # Require at least 5 amide bonds
         return False, f"Found {len(amide_matches)} amide bonds, need at least 5"

    # Check for alpha carbon, stricter definition, must be chiral
    alpha_carbon_pattern = Chem.MolFromSmarts("[C@H]([N])[C](=[O])[C]")
    alpha_carbon_matches = mol.GetSubstructMatches(alpha_carbon_pattern)
    if len(alpha_carbon_matches) < 3:  # Must have at least 3 aminoacids
        return False, f"Found {len(alpha_carbon_matches)} alpha carbon patterns, need at least 3"


    # Check for peptide backbone pattern
    backbone_pattern = Chem.MolFromSmarts("[CX3](=O)[NX2,NX1]-[CX4,CX3]")
    backbone_matches = mol.GetSubstructMatches(backbone_pattern)
    if len(backbone_matches) < 3:
        return False, f"Found {len(backbone_matches)} backbone patterns, need at least 3"

    # Check for a reasonable number of rotatable bonds to see if it can be a peptide
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 8 :
         return False, "Too few rotatable bonds for peptide antibiotic"
    
    # Check for molecular weight (peptide antibiotics are usually > 500)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for peptide antibiotic"
    
    # Check for typical amino acid pattern
    amino_acid_pattern = Chem.MolFromSmarts("[N][C@H](C(=O))[C]") #Check for chiral C
    amino_acid_matches = mol.GetSubstructMatches(amino_acid_pattern)
    if len(amino_acid_matches) <3:
      return False, f"Found {len(amino_acid_matches)} amino acids, need at least 3"
    

    return True, "Contains multiple amide bonds and amino acid structure with a relatively high number of rotatable bonds"