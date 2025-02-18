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

    # Define SMARTS pattern for alpha-amino acid ester structure
    # Pattern: Nitrogen connected to a carbon (alpha), which is connected to carbonyl of ester
    # [NX3] is any nitrogen with three bonds (including implicit Hs)
    # [CX4] is the alpha carbon with four bonds (single bonds to N, C=O, and two others)
    # C(=O)O is the ester carbonyl and oxygen
    amino_ester_pattern = Chem.MolFromSmarts("[NX3]-[CX4]-C(=O)-[OX2]-[#6]")
    
    if mol.HasSubstructMatch(amino_ester_pattern):
        return True, "Contains amino group adjacent to ester carbonyl (alpha position)"
    
    # Check for alternative pattern where alpha carbon is part of a ring (e.g., proline derivatives)
    # This accounts for secondary amines where N is in a ring structure
    ring_amino_ester_pattern = Chem.MolFromSmarts("[NX3&R]-[CX4&R]-C(=O)-[OX2]-[#6]")
    if mol.HasSubstructMatch(ring_amino_ester_pattern):
        return True, "Contains ring amino group adjacent to ester carbonyl (alpha position)"
    
    # If neither pattern matches
    return False, "No alpha-amino acid ester substructure found"