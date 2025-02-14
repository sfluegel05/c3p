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

    # Relaxed peptide bond pattern:
    peptide_bond_pattern = Chem.MolFromSmarts("[CX3](=[OX1])N[X]")
    
    if not mol.HasSubstructMatch(peptide_bond_pattern):
        return False, "Molecule does not contain a peptide bond"

    # Amino acid substructure
    amino_acid_pattern = Chem.MolFromSmarts("[NX3;H2,H1;!$(N*~[C,N,O,S])]C([CX4])(=O)[OX2]")

    amino_acid_matches = mol.GetSubstructMatches(amino_acid_pattern)

    if len(amino_acid_matches) != 2:
         return False, f"Molecule contains {len(amino_acid_matches)} amino-acid residues, should be exactly 2"

    return True, "Contains two amino acid residues connected by a single peptide bond"