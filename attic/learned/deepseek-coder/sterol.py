"""
Classifies: CHEBI:15889 sterol
"""
"""
Classifies: CHEBI:15889 sterol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_sterol(smiles: str):
    """
    Determines if a molecule is a sterol based on its SMILES string.
    A sterol is a 3-hydroxy steroid with a side chain at the 17-position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sterol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the steroid nucleus (four fused rings)
    steroid_nucleus_pattern = Chem.MolFromSmarts("[C@]12[C@]3([C@]4([C@](CC1)(CC[C@@H]4[C@@H]3CC2)[H])[H])[H]")
    if not mol.HasSubstructMatch(steroid_nucleus_pattern):
        return False, "No steroid nucleus found"

    # Check for a hydroxyl group at the 3-position
    hydroxyl_pattern = Chem.MolFromSmarts("[C@]12[C@]3([C@]4([C@](CC1)(CC[C@@H]4[C@@H]3CC2)[H])[H])[H].[OH]")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No hydroxyl group at the 3-position"

    # Check for a side chain at the 17-position
    side_chain_pattern = Chem.MolFromSmarts("[C@]12[C@]3([C@]4([C@](CC1)(CC[C@@H]4[C@@H]3CC2)[H])[H])[H].[C@H](C)CCCC(C)C")
    if not mol.HasSubstructMatch(side_chain_pattern):
        return False, "No side chain at the 17-position"

    # Additional checks for sterol-like properties
    # Check for at least 27 carbons (typical for sterols)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 27:
        return False, "Too few carbons for a sterol"

    # Check for at least 1 oxygen (hydroxyl group)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 1:
        return False, "No oxygen atoms found"

    return True, "Contains a steroid nucleus with a hydroxyl group at the 3-position and a side chain at the 17-position"