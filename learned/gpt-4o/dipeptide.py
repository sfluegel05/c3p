"""
Classifies: CHEBI:46761 dipeptide
"""
from rdkit import Chem

def is_dipeptide(smiles: str):
    """
    Determines if a molecule is a dipeptide based on its SMILES string.
    A dipeptide consists of two amino acid residues linked by a peptide bond.

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
    
    # SMARTS patterns for peptide bond and amino acid residue features
    peptide_bond_pattern = Chem.MolFromSmarts("C(=O)N")
    amine_group_pattern = Chem.MolFromSmarts("[NX3][CX4]")  # e.g., primary amine
    carboxyl_group_pattern = Chem.MolFromSmarts("C(=O)O")
    
    # Look for the peptide bond pattern
    if mol.HasSubstructMatch(peptide_bond_pattern):
        # Locate the peptide bond(s)
        peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)
        
        # Check around each peptide bond for amino acid signatures
        amino_acid_count = 0
        for match in peptide_bond_matches:
            # Ensure there's an amino acid before and after the bond
            if mol.HasSubstructMatch(amine_group_pattern) and mol.HasSubstructMatch(carboxyl_group_pattern):
                amino_acid_count += 1

        # A dipeptide should have exactly two amino acid moieties connected via peptide bond(s)
        if amino_acid_count == 2:
            return True, "Consists of two amino acids connected by a peptide bond"
        else:
            return False, f"Expected 2 amino acid residues, found {amino_acid_count} connections"
            
    else:
        return False, "No peptide bond found or not connected correctly"

# Example usage
# smile = "CC[C@H](C)[C@H](NC(=O)[C@@H](N)CCSC)C(O)=O"  # Example SMILES of Met-Ile
# result = is_dipeptide(smile)
# print(result)