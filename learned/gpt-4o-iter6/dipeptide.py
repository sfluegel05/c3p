"""
Classifies: CHEBI:46761 dipeptide
"""
from rdkit import Chem

def is_dipeptide(smiles: str):
    """
    Determines if a molecule is a dipeptide based on its SMILES string.
    A dipeptide is defined as two amino-acid residues connected by peptide linkages.

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
    
    # Pattern for peptide bond
    peptide_bond_pattern = Chem.MolFromSmarts("N[C;D3](=O)")  # Assuming a peptide bond has an amide N attached to a carbonyl carbon
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    
    # Check if exactly one peptide bond is present
    if len(peptide_bond_matches) != 1:
        return False, f"Expected exactly one peptide bond, found {len(peptide_bond_matches)}"
    
    # Pattern for an amino acid residue (N and C(=O) with an alpha carbon)
    amino_acid_pattern = Chem.MolFromSmarts("N[C;D3][C;D3](=O)") 
    amino_acid_matches = mol.GetSubstructMatches(amino_acid_pattern)
    
    # Check for two amino acid residues
    if len(amino_acid_matches) != 2:
        return False, f"Expected two amino acid residues, found {len(amino_acid_matches)}"

    # Check connectivity of identified patterns
    atoms_match = set()
    for match in peptide_bond_matches + amino_acid_matches:
        atoms_match.update(match)

    # Ensure we aren't counting randomized matches; check for expected dipeptide connectivity
    if len(atoms_match) > len(peptide_bond_matches[0]) + 2 * len(amino_acid_matches[0]) - 1:
        return False, "Patterns found do not form a connected dipeptide with correct linkage"

    return True, "Contains two amino acid residues connected by one peptide bond"

# Test examples
example_smiles = "C[C@@H](O)[C@H](N)C(=O)N[C@@H](Cc1cnc[nH]1)C(O)=O"  # Example: Thr-His
result, reason = is_dipeptide(example_smiles)
print(f"Result: {result}, Reason: {reason}")