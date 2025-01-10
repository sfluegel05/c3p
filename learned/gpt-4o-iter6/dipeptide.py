"""
Classifies: CHEBI:46761 dipeptide
"""
from rdkit import Chem

def is_dipeptide(smiles: str):
    """
    Determines if a molecule is a dipeptide based on its SMILES string.
    A dipeptide contains two amino-acid residues connected by one peptide bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dipeptide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Pattern for peptide bond (CONH linkage)
    peptide_bond_pattern = Chem.MolFromSmarts("N[C;D3](=O)[CX3]")
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    
    # Check if exactly one peptide bond is present
    if len(peptide_bond_matches) != 1:
        return False, f"Expected exactly one peptide bond, found {len(peptide_bond_matches)}"

    # Pattern for an amino acid (NH-CHR-C(=O))
    amino_acid_pattern = Chem.MolFromSmarts("N[C;D3][C;D3](=O)O")
    amino_acid_matches = mol.GetSubstructMatches(amino_acid_pattern)
    
    # Check for exactly two amino acid residues
    if len(amino_acid_matches) < 2:
        return False, f"Expected at least two amino acid residues, found {len(amino_acid_matches)}"

    # Check connectivity of identified patterns
    # Ensure all identified parts form a single connected pattern
    bond_atoms = set(sum(peptide_bond_matches, ()))
    residue_atoms = set(sum(amino_acid_matches, ()))
    
    if len(bond_atoms.intersection(residue_atoms)) < 2:
        return False, "Patterns found do not form a connected dipeptide with correct linkage"

    return True, "Contains two amino acid residues connected by one peptide bond"

# Test examples
example_smiles = "C[C@@H](O)[C@H](N)C(=O)N[C@@H](Cc1cnc[nH]1)C(O)=O"  # Example: Thr-His
result, reason = is_dipeptide(example_smiles)
print(f"Result: {result}, Reason: {reason}")