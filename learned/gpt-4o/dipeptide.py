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

    # SMARTS pattern to detect amino acid residue (generic case)
    alpha_amino_acid_pattern = Chem.MolFromSmarts("[NX3][CX4](C)[CX3](=O)[O,X1,X2]")  # Flexibility for terminal variations

    # Match pattern in molecule
    amino_acid_matches = mol.GetSubstructMatches(alpha_amino_acid_pattern)
    
    # Check if there are exactly two likely alpha amino acid residues
    if len(amino_acid_matches) >= 2:
        # Check for peptide bonds between these residues
        peptide_bond_pattern = Chem.MolFromSmarts("[NX3][CX3](=O)")
        peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)
        
        # Ensure there's at least one peptide bond
        if len(peptide_bond_matches) >= 1:
            return True, "Consists of two amino acids connected by a peptide bond"
        else:
            return False, f"Unexpected number of peptide bonds: {len(peptide_bond_matches)}"
    
    else:
        return False, f"Expected at least 2 amino acid residues, found {len(amino_acid_matches)}"

# Example usage
# smile = "CC[C@H](C)[C@H](NC(=O)[C@@H](N)CCSC)C(O)=O"  # Example SMILES of Met-Ile
# result = is_dipeptide(smile)
# print(result)