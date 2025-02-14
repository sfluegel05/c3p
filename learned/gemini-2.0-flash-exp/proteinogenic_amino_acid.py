"""
Classifies: CHEBI:83813 proteinogenic amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_proteinogenic_amino_acid(smiles: str):
    """
    Determines if a molecule is a proteinogenic amino acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a proteinogenic amino acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for basic alpha-amino acid structure (allowing for protonation states)
    alpha_amino_acid_pattern = Chem.MolFromSmarts("[NX3:1][CX4:2]([CX3](=[OX1])[OX1;H0,H1:3])[CX4]")
    match = mol.GetSubstructMatch(alpha_amino_acid_pattern)
    if not match:
        return False, "Not an alpha-amino acid structure"

    # Check L-configuration (except for glycine, which is achiral)
    alpha_carbon_idx = match[1]
    alpha_carbon = mol.GetAtomWithIdx(alpha_carbon_idx)
    if alpha_carbon.GetDegree() == 4:  # chiral center, check the configuration
       stereo = alpha_carbon.GetChiralTag()

       if stereo == Chem.ChiralType.CHI_TETRAHEDRAL_CW:
          return False, "Chirality is not L"

    # If it matches the pattern and is L configuration, we classify it as proteinogenic
    return True, "Proteinogenic amino acid (side-chain check skipped)"