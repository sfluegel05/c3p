"""
Classifies: CHEBI:15841 polypeptide
"""
"""
Classifies: CHEBI:25004 polypeptide
A peptide containing ten or more amino acid residues.
"""
from rdkit import Chem
from rdkit.Chem import MolFromSmiles, MolFromSmarts
from rdkit.Chem.rdchem import Mol

def is_polypeptide(smiles: str):
    """
    Determines if a molecule is a polypeptide based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polypeptide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define peptide backbone SMARTS pattern
    peptide_backbone_pattern = MolFromSmarts("[N;X3;H2,H1]-,=&@[C;X3](-,=&@[N;X3;H2,H1])-,=&@[C;X3](=O)-,=&@[N;X3;H1]")

    # Check for peptide backbone
    has_peptide_backbone = mol.HasSubstructMatch(peptide_backbone_pattern)
    if not has_peptide_backbone:
        return False, "No peptide backbone detected"

    # Count amino acid residues
    amino_acid_pattern = MolFromSmarts("[N;X3;H2]-,=&@[C;X3](-,=&@[N;X3;H1])-,=&@[C;X3](=O)-,=&@[N;X3;H1]")
    num_residues = len(mol.GetSubstructMatches(amino_acid_pattern))
    if num_residues < 10:
        return False, f"Only {num_residues} amino acid residues, need at least 10"

    return True, f"Contains {num_residues} amino acid residues in a peptide backbone"