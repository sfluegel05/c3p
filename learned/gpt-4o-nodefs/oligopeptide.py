"""
Classifies: CHEBI:25676 oligopeptide
"""
from rdkit import Chem

def is_oligopeptide(smiles: str):
    """
    Determines if a molecule is an oligopeptide based on its SMILES string.
    An oligopeptide can contain peptides with few amino acid residues, starting from dipeptides.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an oligopeptide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES to molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a peptide bond pattern: Carbonyl group adjacent to a nitrogen (amide linkage)
    peptide_bond_pattern = Chem.MolFromSmarts("C(=O)N")

    # Find all instances of peptide bond patterns (-CONH-)
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)

    # Adjust criterion: consider starting with 2 peptide bonds as an oligopeptide
    if len(peptide_bond_matches) >= 2:
        return True, f"Contains {len(peptide_bond_matches)} peptide bonds, classified as oligopeptide"
    
    return False, f"Contains {len(peptide_bond_matches)} peptide bonds, not enough for oligopeptide classification"

# Test the function with an example peptide SMILES
example_smiles = "N[C@@H](C(=O)N[C@@H](C(=O)O)CS)CC1=CC=CC=C1"
result, reason = is_oligopeptide(example_smiles)
print(f"Result: {result}, Reason: {reason}")