"""
Classifies: CHEBI:83813 proteinogenic amino acid
"""
from rdkit import Chem

def is_proteinogenic_amino_acid(smiles: str):
    """
    Determines if a molecule is a proteinogenic amino acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a proteinogenic amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the alpha-amino acid pattern: C-C(N)-C(=O)O
    aa_pattern = Chem.MolFromSmarts("[C@@H](N)C(=O)O")  # Chiral center included
    aa_glycine_pattern = Chem.MolFromSmarts("NC(=O)O")  # Non-chiral for glycine
    
    # Check for chirality 
    has_aa_pattern = mol.HasSubstructMatch(aa_pattern)
    is_glycine = mol.HasSubstructMatch(aa_glycine_pattern)
    
    if not (has_aa_pattern or is_glycine):
        return False, "No proteinogenic amino acid pattern found"
    
    # Check if it has a chiral center labeled with L configuration (natural form)
    if has_aa_pattern:
        for atom in mol.GetAtoms():
            if atom.GetChiralTag() in (Chem.CHI_TETRAHEDRAL_CCW, Chem.CHI_TETRAHEDRAL_CW):
                # Confirm that chiral configuration is L
                chiral_tag = atom.GetChiralTag()
                if chiral_tag != Chem.CHI_TETRAHEDRAL_CCW:
                    return False, "Non-L configuration found for chiral center"
    
    return True, "Pattern matches a proteinogenic amino acid"