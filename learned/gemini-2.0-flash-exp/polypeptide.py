"""
Classifies: CHEBI:15841 polypeptide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_polypeptide(smiles: str):
    """
    Determines if a molecule is a polypeptide based on its SMILES string.
    A polypeptide is defined as a peptide containing ten or more amino acid residues.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polypeptide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS pattern to find a peptide bond
    peptide_bond_pattern = Chem.MolFromSmarts("[-N]-C(=O)-")
    peptide_bonds = mol.GetSubstructMatches(peptide_bond_pattern)
    
    # SMARTS pattern to find amino acid residues within peptide chains
    # This pattern focuses on the N-C-C=O backbone and the alpha carbon with a substitution
    amino_acid_pattern = Chem.MolFromSmarts("[N;H1]-C([*])(-[CX4])C(=O)-")
    amino_acid_matches = mol.GetSubstructMatches(amino_acid_pattern)

    num_amino_acids = len(amino_acid_matches)
    num_peptide_bonds = len(peptide_bonds)
    
    # Check for the minimum number of amino acids.
    if num_amino_acids < 10:
        return False, f"Molecule contains {num_amino_acids} amino acid residues, which is less than the required 10."
    
    # Check that the number of peptide bonds is appropriate for a polypeptide.
    # In most cases, the number of peptide bonds is the number of amino acids minus 1.
    # We allow a tolerance here for cyclic peptides, etc.
    if num_peptide_bonds < num_amino_acids -2:
        return False, f"Molecule contains {num_peptide_bonds} peptide bonds, which seems too low for {num_amino_acids} residues."

    return True, f"Molecule contains {num_amino_acids} amino acid residues connected by peptide bonds, consistent with a polypeptide."