"""
Classifies: CHEBI:83813 proteinogenic amino acid
"""
from rdkit import Chem

def is_proteinogenic_amino_acid(smiles: str):
    """
    Determines if a molecule is a proteinogenic amino acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is a proteinogenic amino acid, False otherwise.
        str: Reason for classification.
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Define a pattern for an alpha-amino acid with allowance for some isotopic labels
    aa_pattern = Chem.MolFromSmarts("[C@@H](N)([CX4,CX3])C(=O)O")
    glycine_pattern = Chem.MolFromSmarts("N[CH2][C](=O)O")  # Generalized glycine pattern
    cys_pattern = Chem.MolFromSmarts("SC[C@H](N)C(=O)O")     # Specific pattern for L-cysteine

    # Check for amino-acid-like basic structure
    has_aa_pattern = mol.HasSubstructMatch(aa_pattern)
    is_glycine = mol.HasSubstructMatch(glycine_pattern)
    is_cysteine = mol.HasSubstructMatch(cys_pattern)
    
    if not (has_aa_pattern or is_glycine or is_cysteine):
        return False, "No proteinogenic amino acid pattern found"

    # Additional length check for chains to avoid peptides or longer derivatives mistakenly classified
    # Acceptable if C count matches common proteinogenic amino acids (generally not exceeding side chain limits for STANDALONE entities)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count > 12:  # Reasonable up to pyrollsine and other larger side chains
        return False, "Too many carbons indicating potential peptide or complex structure"

    # Check chirality, recognizing isotopic chiral configurations when present
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    correct_chirality = False
    has_correct_chirality = False
    for idx, center in chiral_centers:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetSymbol() == 'C':
            if (center == 'S' and not is_cysteine) or (is_cysteine and center == 'R') or is_glycine:
                has_correct_chirality = True
    
    if has_correct_chirality or is_glycine:
        return True, "Pattern matches proteinogenic amino acid"
    
    return False, "Non-L configuration found for chiral center or unknown modifications"