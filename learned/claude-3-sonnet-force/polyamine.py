"""
Classifies: CHEBI:88061 polyamine
"""
"""
Classifies: CHEBI:35586 polyamine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_polyamine(smiles: str):
    """
    Determines if a molecule is a polyamine based on its SMILES string.
    A polyamine is an organic amino compound containing two or more amino groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyamine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count the number of amino groups (-NH2, -NHR, -NR2)
    amino_pattern_1 = Chem.MolFromSmarts("[N;H2&!H0][C;!$(C=[O,N,P,S])]")
    amino_pattern_2 = Chem.MolFromSmarts("[N;H1&!H0;!$(NC=[O,N,P,S])]")
    amino_pattern_3 = Chem.MolFromSmarts("[N;H0&!H1;!$(N([C;!$(C=[O,N,P,S])])(C))]")
    amino_count = len(mol.GetSubstructMatches(amino_pattern_1)) + \
                  len(mol.GetSubstructMatches(amino_pattern_2)) + \
                  len(mol.GetSubstructMatches(amino_pattern_3))

    # Check for at least two amino groups
    if amino_count >= 2:
        # Exclude molecules containing amide groups
        amide_pattern = Chem.MolFromSmarts("C(=O)N")
        if mol.HasSubstructMatch(amide_pattern):
            return False, f"Contains amide group(s), not a polyamine"
        else:
            return True, f"Found {amino_count} amino groups"
    else:
        return False, f"Found only {amino_count} amino group(s), need at least 2"