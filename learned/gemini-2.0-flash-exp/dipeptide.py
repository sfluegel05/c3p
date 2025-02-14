"""
Classifies: CHEBI:46761 dipeptide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_dipeptide(smiles: str):
    """
    Determines if a molecule is a dipeptide based on its SMILES string.
    A dipeptide is a molecule consisting of two amino acid residues connected by a peptide bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dipeptide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Peptide bond pattern
    peptide_bond_pattern = Chem.MolFromSmarts("[CX3](=[OX1])N")
    peptide_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    if len(peptide_matches) < 1:
        return False, "No peptide bond found"


    #Amino acid residue pattern
    amino_acid_pattern = Chem.MolFromSmarts("[NX3][CX4]([HX1])[CX3](=[OX1])[OX2]")
    amino_acid_matches = mol.GetSubstructMatches(amino_acid_pattern)
    if len(amino_acid_matches) < 2:
        return False, f"Found {len(amino_acid_matches)} amino acid residue(s), need at least 2"

    # Check if each aminoacid is connected by a peptide bond:
    peptide_connection = False
    for match in amino_acid_matches:
      for bond in mol.GetBonds():
        if bond.GetBeginAtomIdx() in match and bond.GetEndAtomIdx() in peptide_matches[0]: # check if any amino acids are adjacent to our peptide bond
          peptide_connection = True
    if not peptide_connection:
        return False, "Amino acids not connected via peptide bond"
    
    return True, "Contains two amino acid residues connected by a peptide bond"