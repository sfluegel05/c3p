"""
Classifies: CHEBI:46874 alpha-amino acid ester
"""
"""
Classifies: CHEBI:131313 alpha-amino acid ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_alpha_amino_acid_ester(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid ester based on its SMILES string.
    An alpha-amino acid ester is the formal condensation product of an alpha-amino acid and an alcohol,
    characterized by an ester group (-COO-) and an amino group (-NH2/-NH-) on the adjacent (alpha) carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-amino acid ester, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for alpha-amino acid ester core structure:
    # [NH2/NH]-C(alpha)-C(=O)-O-R
    # Where:
    # - N is primary/secondary amine (not amide)
    # - C(alpha) is connected to both N and ester carbonyl
    # - Ester oxygen connected to any carbon (R group)
    amino_ester_smart = Chem.MolFromSmarts(
        "[NX3;H2,H1;!$(N-C=O)]-!@[CX4H0,CX4H1]"  # Amino group (NH2/NH) connected to alpha carbon
        "(-[#6])"  # Alpha carbon must have at least one substituent (like amino acid side chain)
        "-C(=O)-[OX2]-[#6]"  # Ester group: carbonyl connected to oxygen and R group
    )

    # Check for matches considering possible stereochemistry
    if mol.HasSubstructMatch(amino_ester_smart):
        return True, "Contains amino group adjacent to ester carbonyl (alpha position)"

    # Special case for cyclic amino acids (e.g., proline derivatives)
    # Allows alpha carbon in a ring system
    cyclic_amino_ester_smart = Chem.MolFromSmarts(
        "[NX3;R;!$(N-C=O)]"  # Ring nitrogen (e.g., proline-like)
        "-@[CX4;R]"  # Alpha carbon in ring connected to N
        "-C(=O)-[OX2]-[#6]"  # Ester group
    )
    if mol.HasSubstructMatch(cyclic_amino_ester_smart):
        return True, "Contains cyclic amino group adjacent to ester carbonyl (alpha position)"

    # Check for ester groups in molecule (avoid false positives where carbonyl isn't ester)
    ester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])-[OX2]-[#6]")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester groups present"

    # Final check with relaxed pattern for complex cases (e.g., branched substituents)
    relaxed_smart = Chem.MolFromSmarts(
        "[NX3;H2,H1;!$(N-C=O)]-*!@[CX4]"  # More flexible alpha carbon connection
        "-C(=O)-[OX2]-[#6]"
    )
    if mol.HasSubstructMatch(relaxed_smart):
        return True, "Contains amino group adjacent to ester carbonyl (alpha position)"

    return False, "No alpha-amino acid ester substructure found"