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
    peptide_linkage_pattern = Chem.MolFromSmarts("N[C;X4][C;X3](=O)N")  # Focus on peptide linkage via nitrogen-alpha carbon pattern

    # SMARTS pattern for amino acids: generic amino acid backbone pattern
    amino_acid_pattern = Chem.MolFromSmarts("[NX3][CX4][CX3](=O)[NX3]")

    # Match peptide bonds
    peptide_matches = mol.GetSubstructMatches(peptide_linkage_pattern)
    num_peptide_bonds = len(peptide_matches)

    # Match amino acids
    amino_acid_matches = mol.GetSubstructMatches(amino_acid_pattern)
    num_amino_acids = len(amino_acid_matches)
    
    # For a tetrapeptide, expect 4 residues linked by exactly 3 peptide bonds
    if num_amino_acids >= 4 and num_peptide_bonds == 3:
        return True, "Contains four amino-acid residues connected by peptide linkages"

    return False, f"Found {num_amino_acids} residues and {num_peptide_bonds} peptide bonds, requires 4 residues and 3 peptide bonds"

# Example usage for debug
smiles = "C[C@H](N)C(=O)N[C@@H](Cc1ccccc1)C(=O)N[C@@H](CC(O)=O)C(=O)N[C@@H](CC(O)=O)C(O)=O"
result = is_tetrapeptide(smiles)
print(result)