"""
Classifies: CHEBI:31488 N-acylsphinganine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_N_acylsphinganine(smiles: str):
    """
    Determines if a molecule is an N-acylsphinganine based on its SMILES string.
    An N-acylsphinganine consists of a sphinganine backbone with a fatty acyl group attached to the nitrogen.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acylsphinganine, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Combined SMARTS pattern for sphinganine backbone AND the N-acyl group
    # Includes explicit chain length for fatty acid of 7 C atoms
    sphinganine_acyl_pattern = Chem.MolFromSmarts("[CH2X4][CHX4](O)[CHX4]([CH2X4]O)N[CX3](=[OX1])-[CX4,CX3](-[CX4,CX3])-[CX4,CX3]-[CX4,CX3]-[CX4,CX3]-[CX4,CX3]-[CX4,CX3]")
    if not mol.HasSubstructMatch(sphinganine_acyl_pattern):
        return False, "No N-acylsphinganine substructure found (sphinganine directly connected to acyl group)"
    
    # check that only one such substructure is present.
    matches = mol.GetSubstructMatches(sphinganine_acyl_pattern)
    if len(matches) != 1:
        return False, f"Found {len(matches)} N-acylsphinganine substructures, need exactly 1."

    # Count Carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 14: # a minimum number of carbons in sphinganine backbone
        return False, "Too few carbons for a N-acylsphinganine"
    
    # Count Oxygens
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 3:
      return False, "Too few oxygens"
    if o_count > 8:
        return False, "Too many oxygens"

    return True, "Contains sphinganine backbone with a fatty acyl group attached to the nitrogen"