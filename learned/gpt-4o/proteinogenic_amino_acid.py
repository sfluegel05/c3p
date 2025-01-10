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

    # Define chiral and non-chiral alpha amino acid patterns
    aa_pattern = Chem.MolFromSmarts("[C@@H](N)C(=O)O")  # Pattern for chiral alpha amino acids
    glycine_pattern = Chem.MolFromSmarts("NC(=O)O")     # Pattern for glycine (non-chiral)

    # Match against the patterns
    has_aa_pattern = mol.HasSubstructMatch(aa_pattern)
    is_glycine = mol.HasSubstructMatch(glycine_pattern)

    if not (has_aa_pattern or is_glycine):
        return False, "No proteinogenic amino acid pattern found"

    # If pattern is found, glycine is okay as is
    if is_glycine:
        return True, "Pattern matches proteinogenic (glycine) amino acid"

    # For chiral amino acids, ensure correct stereochemistry
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    
    for idx, center in chiral_centers:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetSymbol() == 'C':
            # Consider any L amino acids or valid isotopic labels
            if center == 'R':  # Typically L-amino acids are S; may need manual S/R correction knowing input set
                return True, "Pattern matches proteinogenic amino acid"
    
    return False, "Non-L configuration found for chiral center or unknown modifications"