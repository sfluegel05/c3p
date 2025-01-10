"""
Classifies: CHEBI:25409 monoterpenoid
"""
"""
Classifies: monoterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monoterpenoid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid based on its SMILES string.
    A monoterpenoid is a terpenoid derived from a monoterpene (C10 skeleton).
    The term includes compounds in which the C10 skeleton of the parent monoterpene
    has been rearranged or modified by the removal of one or more skeletal atoms
    (generally methyl groups).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count the number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Monoterpenoids typically have 7 to 10 carbons (due to possible removal of methyl groups)
    if c_count < 7 or c_count > 10:
        return False, f"Number of carbons ({c_count}) not in expected range for monoterpenoids (7-10)"
    
    # Attempt to identify isoprene units (C5H8)
    # Define isoprene substructure pattern
    isoprene_smarts = "[CH2]=[CH][CH]=[CH2]"  # Simplified pattern for isoprene unit
    isoprene_pattern = Chem.MolFromSmarts(isoprene_smarts)
    if not mol.HasSubstructMatch(isoprene_pattern):
        # Check for branched isoprene units
        branched_isoprene_smarts = "C(=C)C=C"
        branched_isoprene_pattern = Chem.MolFromSmarts(branched_isoprene_smarts)
        if not mol.HasSubstructMatch(branched_isoprene_pattern):
            return False, "No isoprene units found"

    # Check for common terpenoid functional groups
    functional_group_patterns = [
        "[OX2H]",        # Hydroxyl group
        "[CX3]=[OX1]",   # Carbonyl group (ketone or aldehyde)
        "[CX3](=O)[OX2H1]",  # Carboxylic acid
        "[CX3](=O)[OX2][CX4]",  # Ester
        "[OX2][CX4]",    # Ether
    ]
    fg_found = False
    for fg_smarts in functional_group_patterns:
        fg_pattern = Chem.MolFromSmarts(fg_smarts)
        if mol.HasSubstructMatch(fg_pattern):
            fg_found = True
            break
    if not fg_found:
        return False, "No typical terpenoid functional groups found"

    return True, "Molecule has characteristics of a monoterpenoid"