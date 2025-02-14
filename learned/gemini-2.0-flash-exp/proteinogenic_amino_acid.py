"""
Classifies: CHEBI:83813 proteinogenic amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_proteinogenic_amino_acid(smiles: str):
    """
    Determines if a molecule is a proteinogenic amino acid based on its SMILES string.
    This version checks for the core alpha-amino acid structure and the L configuration but skips side-chain specific checks.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a proteinogenic amino acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more specific SMARTS pattern for the alpha-amino acid core, allowing for protonation.
    # The carboxyl oxygen can be single bonded (O-) or double bonded (=O) and with or without H atoms.
    # The nitrogen can be protonated or not.
    alpha_amino_acid_pattern = Chem.MolFromSmarts("[NX3;H0,H1,H2,H3:1][CX4:2]([CX3:3](=[OX1])[OX1;H0,H1:4])[CX4]")
    match = mol.GetSubstructMatch(alpha_amino_acid_pattern)

    if not match:
        return False, "Not an alpha-amino acid structure"


    # Check L-configuration (except for glycine, which is achiral).
    # Glycine's alpha carbon is a CH2, so it will have only 3 neighbors.
    alpha_carbon_idx = match[1]
    alpha_carbon = mol.GetAtomWithIdx(alpha_carbon_idx)

    # Only check chirality if alpha carbon has 4 neighbors.
    if alpha_carbon.GetDegree() == 4:
        stereo = alpha_carbon.GetChiralTag()
        if stereo == Chem.ChiralType.CHI_TETRAHEDRAL_CW: #L-configuration is CCW
            return False, "Chirality is not L"
            

    return True, "Proteinogenic amino acid (side-chain check skipped)"