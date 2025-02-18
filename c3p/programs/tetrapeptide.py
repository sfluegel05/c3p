"""
Classifies: CHEBI:48030 tetrapeptide
"""
from rdkit import Chem

def is_tetrapeptide(smiles: str):
    """
    Determines if a molecule is a tetrapeptide based on its SMILES string.
    A tetrapeptide is any molecule that contains four amino-acid residues connected by peptide linkages.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrapeptide, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for a peptide linkage
    peptide_linkage_pattern = Chem.MolFromSmarts("[NX3][CX3](=O)[NX3X2]")  # N-any C(=O)-N pattern for amide linkages

    # Pattern for standard amino acid structure
    amino_acid_patterns = [
        Chem.MolFromSmarts("[NX3][CX4;H2][CX3](=O)[NX3]"),
        Chem.MolFromSmarts("[NX3][CX4;H][CX3](=O)[NX3]"),
        Chem.MolFromSmarts("[NX3][CX4][CX3](=O)[NX3]"),
    ]
    
    # Match peptide bonds
    peptide_matches = mol.GetSubstructMatches(peptide_linkage_pattern)
    num_peptide_bonds = len(peptide_matches)

    # Match amino acids using all patterns
    num_amino_acids = 0
    for pattern in amino_acid_patterns:
        matches = mol.GetSubstructMatches(pattern)
        num_amino_acids += len(matches)

    # For a tetrapeptide, expect 4 residues linked by exactly 3 peptide bonds
    if num_amino_acids >= 4 and num_peptide_bonds == 3:
        return True, f"Contains {num_amino_acids} amino-acid residues connected by 3 peptide linkages"

    return False, f"Found {num_amino_acids} residues and {num_peptide_bonds} peptide bonds, requires 4 residues and 3 peptide bonds"

# Example usage for debug
smiles = "C[C@@H](O)[C@H](NC(=O)[C@H](C)N)C(=O)N[C@@H](C)C(=O)N1CCC[C@H]1C(O)=O"
result = is_tetrapeptide(smiles)
print(result)