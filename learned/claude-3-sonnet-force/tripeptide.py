"""
Classifies: CHEBI:47923 tripeptide
"""
"""
Classifies: CHEBI:33938 tripeptide
A tripeptide is any oligopeptide that consists of three amino-acid residues connected by peptide linkages.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_tripeptide(smiles: str):
    """
    Determines if a molecule is a tripeptide based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tripeptide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for 3 connected amino acid residues (N-C(=O)-C)
    amino_acid_pattern = Chem.MolFromSmarts("[N&x3]-[C&x3](=O)-[C&x3]")
    amino_acid_matches = mol.GetSubstructMatches(amino_acid_pattern)
    if len(amino_acid_matches) != 3:
        return False, f"Found {len(amino_acid_matches)} amino acid residues, need exactly 3"

    # Check for connected amino acid residues
    connected_residues = []
    for match in amino_acid_matches:
        residue = tuple(match)
        if any(residue[0] in prev_residue and residue[2] in next_residue for prev_residue, next_residue in zip(connected_residues, connected_residues[1:])):
            connected_residues.append(residue)
        else:
            return False, "Amino acid residues not connected by peptide bonds"

    # Check for N-terminus and C-terminus
    n_terminus_pattern = Chem.MolFromSmarts("[N&x3]-[C&x3](=O)-")
    c_terminus_pattern = Chem.MolFromSmarts("-[C&x3](=O)-[O&x1,N&x3]")
    if not mol.HasSubstructMatch(n_terminus_pattern) or not mol.HasSubstructMatch(c_terminus_pattern):
        return False, "Missing N-terminus or C-terminus"

    # Check for common modifications or substituents
    cyclic_pattern = Chem.MolFromSmarts("[N&x3]-[C&x3](=O)-[C&x3]-1-[C&x3]-[C&x3]-[C&x3]-1")
    if mol.HasSubstructMatch(cyclic_pattern):
        reason = "Contains 3 amino acid residues connected by peptide bonds, with a cyclic residue (diketopiperazine)"
    else:
        reason = "Contains 3 amino acid residues connected by peptide bonds"

    return True, reason