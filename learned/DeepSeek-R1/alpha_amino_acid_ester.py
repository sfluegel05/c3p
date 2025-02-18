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

    # Improved SMARTS patterns to account for varied bonding and ring systems
    # Pattern 1: General case (linear or non-ring) with any bond type between N and alpha carbon
    # [NX3] includes primary, secondary amines and aromatic nitrogens with three bonds
    # [CX3,CX4] allows alpha carbon with 3 or 4 bonds (including ring carbons)
    amino_ester_pattern = Chem.MolFromSmarts("[NX3]~[CX3,CX4]-C(=O)-[OX2]-[#6]")
    
    # Pattern 2: Nitrogen in a ring system (e.g., proline derivatives)
    # Checks for ring-context amino groups connected to ester carbonyl
    ring_amino_ester_pattern = Chem.MolFromSmarts("[NX3&R]~[CX3,CX4&R]-C(=O)-[OX2]-[#6]")
    
    # Check both patterns
    if mol.HasSubstructMatch(amino_ester_pattern) or mol.HasSubstructMatch(ring_amino_ester_pattern):
        return True, "Contains amino group adjacent to ester carbonyl (alpha position)"
    
    # Additional check for cases where ester oxygen is part of a larger structure
    # and alpha carbon has substituents (like in branched chains)
    # This pattern allows up to two arbitrary atoms between N and ester carbonyl
    extended_pattern = Chem.MolFromSmarts("[NX3]~[CX3,CX4](-[#6])-C(=O)-[OX2]-[#6]")
    if mol.HasSubstructMatch(extended_pattern):
        return True, "Contains amino group adjacent to ester carbonyl (alpha position)"
    
    return False, "No alpha-amino acid ester substructure found"