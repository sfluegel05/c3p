"""
Classifies: CHEBI:46761 dipeptide
"""
"""
Classifies: CHEBI:46761 dipeptide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_dipeptide(smiles: str):
    """
    Determines if a molecule is a dipeptide based on its SMILES string.
    A dipeptide contains two amino-acid residues connected by a peptide bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dipeptide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find all amide bonds (C(=O)-N with single bond to N)
    amide_pattern = Chem.MolFromSmarts("[CX3](=O)[NX3H1]")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if len(amide_matches) != 1:
        return False, f"Found {len(amide_matches)} amide bonds, need exactly 1"

    # Check that the amide connects two alpha carbons (each with possible side chains)
    c_atom = amide_matches[0][0]
    n_atom = amide_matches[0][1]

    # Get neighboring carbons to the amide C and N (alpha carbons)
    alpha1 = None
    for neighbor in mol.GetAtomWithIdx(c_atom).GetNeighbors():
        if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() != n_atom:
            alpha1 = neighbor
            break
    if not alpha1:
        return False, "No alpha carbon adjacent to amide C"

    alpha2 = None
    for neighbor in mol.GetAtomWithIdx(n_atom).GetNeighbors():
        if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() != c_atom:
            alpha2 = neighbor
            break
    if not alpha2:
        return False, "No alpha carbon adjacent to amide N"

    # Check if each alpha carbon has at least one substituent (side chain)
    # Alpha1 should have at least two bonds (to amide C and another group)
    if len(alpha1.GetBonds()) < 2:
        return False, "First alpha carbon lacks substituents"

    # Alpha2 should have at least two bonds (to amide N and another group)
    if len(alpha2.GetBonds()) < 2:
        return False, "Second alpha carbon lacks substituents"

    return True, "Two amino acid residues connected by a peptide bond"