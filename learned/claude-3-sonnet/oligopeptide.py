"""
Classifies: CHEBI:25676 oligopeptide
"""
"""
Classifies: CHEBI:16670 oligopeptide
"""

from rdkit import Chem
from rdkit.Chem import Descriptors

AA_CODES = "Ala Arg Asn Asp Cys Gln Glu Gly His Ile Leu Lys Met Phe Pro Ser Thr Trp Tyr Val".split()

def is_oligopeptide(smiles: str):
    """
    Determines if a molecule is an oligopeptide based on its SMILES string.
    An oligopeptide is a peptide containing a relatively small number of amino acids.

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

    # Get molecular weight
    mol_wt = Descriptors.ExactMolWt(mol)

    # Count amino acid residues
    long_name = Descriptors.GetMULongName(mol)
    aa_residues = sum(code in long_name for code in AA_CODES)

    # Look for peptide bonds
    peptide_bond_pattern = Chem.MolFromSmarts("[N;X3][C;X3](=[O;X1])[C;X3][N;X3]")
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)

    # Check for common functional groups
    has_carboxyl = mol.HasSubstructMatch(Chem.MolFromSmarts("[C;X3](=O)[O;X1]"))
    has_amino = mol.HasSubstructMatch(Chem.MolFromSmarts("[N;X3;H2,H1]"))
    has_hydroxyl = mol.HasSubstructMatch(Chem.MolFromSmarts("[O;H1]"))
    has_aromatic = mol.HasSubstructMatch(Chem.MolFromSmarts("[a]"))

    # Classify as oligopeptide
    if aa_residues >= 2 and len(peptide_bond_matches) >= 1 and has_carboxyl and has_amino and (has_hydroxyl or has_aromatic):
        if mol_wt > 500 and mol_wt < 5000:
            return True, f"Contains {aa_residues} amino acid residues connected by peptide bonds, with common functional groups"
        else:
            return False, f"Molecular weight ({mol_wt:.2f} Da) outside typical oligopeptide range"
    else:
        return False, "Lacks characteristic features of oligopeptides"