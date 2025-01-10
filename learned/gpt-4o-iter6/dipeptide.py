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
        bool: True if the molecule is a dipeptide, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for two amino acids connected by one peptide bond
    # Pattern: Two segments, each with an amino acid (NH-CH-CO) connected by a peptide bond (CO-NH)
    dipeptide_pattern = Chem.MolFromSmarts("N[C;D3][CX3](=O)N[C;D3][CX3](=O)[O,N]")  # Sequential connection
    
    # Check for the dipeptide pattern in the molecule
    if mol.HasSubstructMatch(dipeptide_pattern):
        return True, "Contains two amino acid residues connected by one peptide bond"
    
    return False, "No connected dipeptide structure found"

# Test examples
example_smiles = "C[C@@H](O)[C@H](N)C(=O)N[C@@H](Cc1cnc[nH]1)C(O)=O"  # Example: Thr-His
result, reason = is_dipeptide(example_smiles)
print(f"Result: {result}, Reason: {reason}")