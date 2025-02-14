"""
Classifies: CHEBI:25676 oligopeptide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_oligopeptide(smiles: str):
    """
    Determines if a molecule is an oligopeptide based on its SMILES string.
    An oligopeptide is a short chain of amino acids linked by peptide bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an oligopeptide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for peptide bonds (-C(=O)N-)
    peptide_bond_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[NX3]")
    peptide_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    num_peptide_bonds = len(peptide_matches)
    if num_peptide_bonds < 1: #An oligopeptide should have at least one peptide bond
        return False, f"Found {num_peptide_bonds} peptide bonds, at least 1 required"

    # Count amino acid residues - must be between 1 and 50
    # Since each peptide bond links two amino acids, the number of residues is simply one more than the number of peptide bonds, however this is an approximation. We will use the number of peptide bonds as a proxy for residues for now.
    num_residues = num_peptide_bonds + 1

    if num_residues < 1 or num_residues > 50:
            return False, f"Found {num_residues} amino acid residues, must be between 1 and 50"

    # Check for at least one carboxyl and at least one amine
    has_terminal_carboxyl = mol.HasSubstructMatch(Chem.MolFromSmarts("C(=O)[OH]")) or mol.HasSubstructMatch(Chem.MolFromSmarts("C(=O)[O-]"))
    has_terminal_amino    = mol.HasSubstructMatch(Chem.MolFromSmarts("[NH2]")) or mol.HasSubstructMatch(Chem.MolFromSmarts("[NH3+]")) or mol.HasSubstructMatch(Chem.MolFromSmarts("[NHX2]"))

    # A molecule should have at least one group that could form a peptide bond
    if not (has_terminal_carboxyl or num_peptide_bonds > 0) or not (has_terminal_amino or num_peptide_bonds > 0):
        return False, "Not a typical oligopeptide structure without terminal amino or carboxyl group and without peptide bonds"
    

    return True, "Has peptide bonds and amino acid residues, within size range of a typical oligopeptide"