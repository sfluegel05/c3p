"""
Classifies: CHEBI:15841 polypeptide
"""
"""
Classifies: CHEBI:25004 polypeptide
A peptide containing ten or more amino acid residues.
"""
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import Fragments

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
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for peptide backbone
    has_peptide_backbone = any(Fragments.fr_peptide.fp_bit in fp for fp in Fragments.get_fragments(mol))
    if not has_peptide_backbone:
        return False, "No peptide backbone detected"

    # Count amino acid residues
    num_residues = Descriptors.NumResidues(mol)
    if num_residues < 10:
        return False, f"Only {num_residues} amino acid residues, need at least 10"

    return True, f"Contains {num_residues} amino acid residues in a peptide backbone"