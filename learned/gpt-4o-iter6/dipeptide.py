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
    
    # Pattern for peptide bond (CO-NH linkage)
    peptide_bond_pattern = Chem.MolFromSmarts("[NX3][CX3](=O)")
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    
    # Check if at least one peptide bond is present
    if len(peptide_bond_matches) < 1:
        return False, f"No peptide bond found in structure"

    # Pattern for an amino acid residues specifically via carboxyl end
    amino_acid_pattern = Chem.MolFromSmarts("N[C;D3][C;D3](=O)O")
    amino_acid_matches = mol.GetSubstructMatches(amino_acid_pattern)
    
    # Initial len check for amino acid residues on carboxyl end
    if len(amino_acid_matches) < 1:
        return False, f"Expected at least one amino acid (carboxyl end), found {len(amino_acid_matches)}"

    # Check for each peptide bond, ensure connection to carboxyl groups
    for bond in peptide_bond_matches:
        for aa_match in amino_acid_matches:
            # Verify connectivity between amino acid carboxyl and peptide bond
            bond_atoms = set(bond)
            residue_atoms = set(aa_match)
            
            # If found a valid connection
            if len(bond_atoms.intersection(residue_atoms)) >= 2:
                return True, "Contains two amino acid residues connected by one peptide bond"

    return False, "Patterns found do not form a connected dipeptide with correct linkage"

# Test examples
example_smiles = "C[C@@H](O)[C@H](N)C(=O)N[C@@H](Cc1cnc[nH]1)C(O)=O"  # Example: Thr-His
result, reason = is_dipeptide(example_smiles)
print(f"Result: {result}, Reason: {reason}")