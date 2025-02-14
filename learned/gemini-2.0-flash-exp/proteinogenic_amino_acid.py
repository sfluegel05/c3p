"""
Classifies: CHEBI:83813 proteinogenic amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_proteinogenic_amino_acid(smiles: str):
    """
    Determines if a molecule is a proteinogenic amino acid based on its SMILES string.
    This version checks for the core alpha-amino acid structure, the L configuration, and for the presence of an appropriate R-group attachment.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a proteinogenic amino acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Updated SMARTS pattern to account for protonated/deprotonated forms and a connecting carbon.
    # This matches the core of any alpha amino acid and includes a carbon atom (or deuterium-substituted carbon) attached to the alpha carbon as the R group connection.
    alpha_amino_acid_pattern = Chem.MolFromSmarts("[NX3;H0,H1,H2,H3:1][CX4:2]([CX3:3](=[OX1;H0,H1:4])[OX1;H0,H1,H2:5])[CX4;D1,D2:6]")
    match = mol.GetSubstructMatch(alpha_amino_acid_pattern)

    if not match:
         return False, "Not an alpha-amino acid structure"
        
    # Check L-configuration (except for glycine, which is achiral)
    alpha_carbon_idx = match[1]
    alpha_carbon = mol.GetAtomWithIdx(alpha_carbon_idx)

    # If alpha carbon has 3 neighbors it's glycine, skip chirality check
    if alpha_carbon.GetDegree() == 4: 
        stereo = alpha_carbon.GetChiralTag()
        if stereo == Chem.ChiralType.CHI_TETRAHEDRAL_CCW: #L-configuration is CCW
           pass #This is correct configuration
        elif stereo == Chem.ChiralType.CHI_TETRAHEDRAL_CW:
           return False, "Chirality is not L"
        elif stereo == Chem.ChiralType.CHI_UNSPECIFIED:
          pass # Chiral center is not specified. Allow this, since e.g. L-proline is allowed
        else:
           return False, "Chirality not defined for alpha carbon"
       
    return True, "Proteinogenic amino acid"