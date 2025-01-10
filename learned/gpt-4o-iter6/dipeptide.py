"""
Classifies: CHEBI:46761 dipeptide
"""
from rdkit import Chem

def is_dipeptide(smiles: str):
    """
    Determines if a molecule is a dipeptide based on its SMILES string.
    A dipeptide contains two amino-acid residues connected by peptide bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a dipeptide, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for two amino acids connected by peptide bonds
    # Generalized pattern considering various configurations and protection groups
    dipeptide_pattern = Chem.MolFromSmarts("[NX3][C;D3][CX3](=O)N[C;D3][CX3](=O)[O,N,C]")  
    
    # Check for the dipeptide pattern in the molecule
    if mol.HasSubstructMatch(dipeptide_pattern):
        return True, "Contains two amino acid residues connected by peptide bonds"

    return False, "No connected dipeptide structure found"

# Example test
example_smiles = "C[C@@H](O)[C@H](N)C(=O)N[C@@H](Cc1cnc[nH]1)C(O)=O"  # Example: Thr-His
result, reason = is_dipeptide(example_smiles)
print(f"Result: {result}, Reason: {reason}")

# Note: The improvements lie in flexibility for terminal and linkage tolerance