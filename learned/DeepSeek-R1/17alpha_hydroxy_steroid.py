"""
Classifies: CHEBI:35342 17alpha-hydroxy steroid
"""
"""
Classifies: CHEBI:17450 17alpha-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_17alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 17alpha-hydroxy steroid based on its SMILES string.
    A 17alpha-hydroxy steroid has a hydroxyl group in the alpha configuration at position 17 of the steroid nucleus.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 17alpha-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the steroid nucleus (four fused rings)
    # This SMARTS pattern approximates the tetracyclic structure of steroids
    steroid_core = Chem.MolFromSmarts("[C@]12[C@@H]3[C@@H]4[C@@H]([C@@H]1CCCC2)[C@H]1CC[C@@H]4[C@@H]31")
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid nucleus detected"
    
    # Check for hydroxyl group at position 17 with alpha configuration
    # SMARTS pattern for 17alpha-OH: hydroxyl attached to a carbon in the D ring with specific stereochemistry
    alpha_oh_pattern = Chem.MolFromSmarts("[C@]([OH])(C)(C)C")
    alpha_matches = mol.GetSubstructMatches(alpha_oh_pattern)
    if not alpha_matches:
        # Check alternative pattern where the hydroxyl is part of a five-membered ring (D ring)
        alpha_oh_ring = Chem.MolFromSmarts("[C@]1(O)CCCC1")
        alpha_ring_matches = mol.GetSubstructMatches(alpha_oh_ring)
        if not alpha_ring_matches:
            return False, "No 17-alpha hydroxyl group found"
    
    return True, "Steroid nucleus with 17-alpha hydroxyl group"