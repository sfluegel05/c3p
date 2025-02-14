"""
Classifies: CHEBI:35346 11beta-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_11beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 11beta-hydroxy steroid based on its SMILES string.
    An 11beta-hydroxy steroid has a steroid core with a beta-configured hydroxyl group at position 11.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an 11beta-hydroxy steroid, False otherwise
        str: Reason for the classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more specific SMARTS for the 11beta-hydroxy group, anchored to part of the steroid core.
    # This pattern looks for a C-C-C chain where the middle C is beta and has a hydroxyl.
    # It also requires that part of the steroid core ring system be present
    hydroxy_11beta_pattern = Chem.MolFromSmarts("[C]([C])([C@H]1[C@@](O)([C])~[C]~[C]1)")


    if not mol.HasSubstructMatch(hydroxy_11beta_pattern):
        return False, "No 11-beta-hydroxy group found in the context of a steroid ring"

    # Define a more flexible SMARTS for the steroid core
    # This SMARTS describes the four fused rings with some flexibility using ~ for any bond
    steroid_core_pattern = Chem.MolFromSmarts("[C,CH,CH2,CH3]~[C,CH,CH2,CH3]~[C,CH,CH2,CH3]~[C,CH,CH2,CH3]~[C,CH,CH2,CH3]1~[C,CH,CH2,CH3]~[C,CH,CH2,CH3]~[C,CH,CH2,CH3]~[C,CH,CH2,CH3]2~[C,CH,CH2,CH3]~[C,CH,CH2,CH3]~[C,CH,CH2,CH3]~[C,CH,CH2,CH3]1~[C,CH,CH2,CH3]~[C,CH,CH2,CH3]2")
    if not mol.HasSubstructMatch(steroid_core_pattern):
         return False, "No steroid core structure detected (no 4 fused rings)"

    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 17:
         return False, f"Too few carbons, expected 17 in core, got {carbon_count}"


    return True, "Molecule has a steroid core with an 11-beta-hydroxyl group"