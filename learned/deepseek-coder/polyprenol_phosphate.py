"""
Classifies: CHEBI:16460 polyprenol phosphate
"""
"""
Classifies: CHEBI:37550 polyprenol phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_polyprenol_phosphate(smiles: str):
    """
    Determines if a molecule is a polyprenol phosphate based on its SMILES string.
    A polyprenol phosphate is a prenol phosphate resulting from the formal condensation
    of the terminal allylic hydroxy group of a polyprenol with 1 mol eq. of phosphoric acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyprenol phosphate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for phosphate group pattern (P with at least one O and one O-)
    phosphate_pattern = Chem.MolFromSmarts("[PX4]([OX2])([OX2])([OX2])([OX1-]?)")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"

    # Look for polyprenol chain (alternating double bonds in a long chain)
    polyprenol_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    polyprenol_matches = mol.GetSubstructMatches(polyprenol_pattern)
    if len(polyprenol_matches) < 1:
        return False, "No polyprenol chain found"

    # Check if the phosphate is attached to the terminal allylic hydroxy group
    # This is a complex pattern, so we'll look for a phosphate attached to a carbon with a double bond
    phosphate_attachment_pattern = Chem.MolFromSmarts("[PX4]([OX2])([OX2])([OX2])([OX1-]?)[CX3]=[CX4]")
    if not mol.HasSubstructMatch(phosphate_attachment_pattern):
        return False, "Phosphate not attached to terminal allylic carbon"

    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Chains too short to be polyprenol"

    # Check molecular weight - polyprenol phosphates typically >400 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, "Molecular weight too low for polyprenol phosphate"

    return True, "Contains polyprenol chain with phosphate group attached to terminal allylic carbon"