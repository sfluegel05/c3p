"""
Classifies: CHEBI:83813 proteinogenic amino acid
"""
from rdkit import Chem
from rdkit.Chem import rdqueries

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
        
    # Define a pattern for an alpha-amino acid - broader to include derivatives or modifications
    aa_pattern = Chem.MolFromSmarts("[C@@H](N)([CX4,CX3])C(=O)O")  # Pattern matching chiral alpha amino acids
    glycine_pattern = Chem.MolFromSmarts("NCC(=O)O")                # Pattern matching glycine backbone
    cys_pattern = Chem.MolFromSmarts("SC[C@@H](N)C(=O)O")           # Specific pattern for L-cysteine

    # Match against the patterns
    has_aa_pattern = mol.HasSubstructMatch(aa_pattern)
    is_glycine = mol.HasSubstructMatch(glycine_pattern)
    is_cysteine = mol.HasSubstructMatch(cys_pattern)
    
    if not (has_aa_pattern or is_glycine or is_cysteine):
        return False, "No proteinogenic amino acid pattern found"

    # Check chiral centers correctly for non-glycine amino acids
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    for idx, center in chiral_centers:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetSymbol() == 'C':
            # Consider L-amino acids must be the S configuration except for cysteine (R)
            # Check correctly set flags or isotopic labels
            if (center == 'S' or is_cysteine and center == 'R') or is_glycine:
                return True, "Pattern matches proteinogenic amino acid"

    return False, "Non-L configuration found for chiral center or unknown modifications"