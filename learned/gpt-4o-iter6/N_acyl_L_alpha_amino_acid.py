"""
Classifies: CHEBI:48927 N-acyl-L-alpha-amino acid
"""
from rdkit import Chem

def is_N_acyl_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acyl-L-alpha-amino acid based on its SMILES string.
    This class is defined by an L-alpha-amino acid backbone with an N-acyl substituent.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acyl-L-alpha-amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (False, "Invalid SMILES string")
    
    # General pattern for L-alpha-amino acid backbone with stereochemistry specified
    l_alpha_amino_acid_pattern = Chem.MolFromSmarts('[C@H](N)(*)C(=O)O')  # Central chiral carbon with amine and carboxyl connected
    if not mol.HasSubstructMatch(l_alpha_amino_acid_pattern):
        return (False, "No chiral L-alpha-amino acid backbone found")

    # Identify Nitrogen having an acyl group directly bonded
    n_acyl_pattern = Chem.MolFromSmarts('N[C](=O)C')  # N attached to carbonyl and chain
    n_acyl_matches = mol.GetSubstructMatches(n_acyl_pattern)

    if not n_acyl_matches:
        return (False, "No N-acyl group found attached to the expected nitrogen")

    # Further confirm attachment of N-acyl to correct nitrogen
    for match in n_acyl_matches:
        # Check if the nitrogen of n_acyl is connected back to the alpha carbon
        nitrogen_idx = match[0]
        BondsWithNGroup = mol.GetAtomWithIdx(nitrogen_idx).GetDegree()  # number of bonds should match amino acid structure
        if BondsWithNGroup == 1:  # Verify isolated N-acyl rather than part of peptide linkage
            return (True, "Contains L-alpha-amino acid backbone with an N-acyl substituent correctly attached")

    return (False, "N-acyl group not correctly integrated into L-alpha-amino acid nitrogen")