"""
Classifies: CHEBI:48039 dihydroflavonols
"""
"""
Classifies: CHEBI:???? dihydroflavonols (3-hydroxyflavanones)
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_dihydroflavonols(smiles: str):
    """
    Determines if a molecule is a dihydroflavonol based on its SMILES string.
    A dihydroflavonol is a hydroxyflavanone with a hydroxyl group at position 3 of the heterocyclic ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dihydroflavonol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for 3-hydroxyflavanone core
    # The pattern matches the flavanone skeleton with a hydroxyl group on position 3
    # Flavanone core is a chroman-4-one (O-C2-C3-OH-C4=O in a six-membered ring)
    # B ring (c2ccccc2) can have any substituents
    pattern = Chem.MolFromSmarts("[OH]C1C(=O)c2ccccc2OC1")
    
    if not mol.HasSubstructMatch(pattern):
        return False, "Does not contain 3-hydroxyflavanone core structure"
    
    # Additional check to exclude false positives: ensure the hydroxyl is directly attached to the C3 in the ring
    # Get the matching atoms to verify connectivity
    matches = mol.GetSubstructMatches(pattern)
    if not matches:
        return False, "No 3-hydroxyflavanone core match"
    
    # Check that the hydroxyl is part of the ring system
    for match in matches:
        hydroxyl_atom = mol.GetAtomWithIdx(match[0])
        ring_info = mol.GetRingInfo()
        if hydroxyl_atom.IsInRing() and any(len(r) >= 6 for r in ring_info.AtomRings()):
            return True, "Contains 3-hydroxyflavanone core structure"
    
    return False, "Hydroxyl not in a six-membered ring (flavanone structure)"