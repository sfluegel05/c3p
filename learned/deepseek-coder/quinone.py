"""
Classifies: CHEBI:36141 quinone
"""
"""
Classifies: CHEBI:36141 quinone
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_quinone(smiles: str):
    """
    Determines if a molecule is a quinone based on its SMILES string.
    A quinone is a compound with a fully conjugated cyclic dione structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a quinone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a cyclic structure
    if not mol.GetRingInfo().NumRings():
        return False, "No cyclic structure found"

    # Look for two carbonyl groups (C=O) in the molecule
    carbonyl_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")
    carbonyl_matches = mol.GetSubstructMatches(carbonyl_pattern)
    if len(carbonyl_matches) < 2:
        return False, f"Found {len(carbonyl_matches)} carbonyl groups, need at least 2"

    # Check if the carbonyl groups are in a 1,4-relationship within the same ring
    quinone_pattern = Chem.MolFromSmarts("[CX3]=[OX1]1[CX3]=[CX3][CX3]=[CX3]1")
    if not mol.HasSubstructMatch(quinone_pattern):
        # Try alternative patterns for polycyclic and heterocyclic quinones
        extended_pattern = Chem.MolFromSmarts("[CX3]=[OX1]1[CX3]=[CX3][CX3]=[CX3][CX3]=[CX3]1")
        if not mol.HasSubstructMatch(extended_pattern):
            return False, "No 1,4-dione structure found"

    # Check for conjugation by ensuring the molecule is planar
    # This is a more reliable indicator of conjugation than simple double bond checks
    if not rdMolDescriptors.CalcNumAromaticRings(mol):
        # If not aromatic, check for extended conjugation
        conjugated_pattern = Chem.MolFromSmarts("[CX3]=[CX3]~[CX3]=[CX3]")
        if not mol.HasSubstructMatch(conjugated_pattern):
            return False, "No conjugated system found"

    return True, "Contains a fully conjugated cyclic dione structure"