"""
Classifies: CHEBI:36836 3beta-hydroxy steroid
"""
"""
Classifies: CHEBI:36804 3beta-hydroxy steroid
"""
from rdkit import Chem

def is_3beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3beta-hydroxy steroid based on its SMILES string.
    A 3beta-hydroxy steroid is a steroid with a hydroxyl group at the 3-position in the beta orientation.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure stereochemistry is assigned
    Chem.AssignAtomChiralTagsFromStructure(mol)
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    
    # Steroid nucleus SMARTS pattern (cyclopentanoperhydrophenanthrene)
    steroid_pattern = Chem.MolFromSmarts('''
        [#6]1[#6][#6][#6]2[#6]([#6]1)[#6][#6][#6]3[#6]2[#6][#6][#6]4[#6]3[#6][#6][#6]([#6]4)
    ''')
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "Steroid core not found"
    
    # 3beta-hydroxy group SMARTS pattern
    # The 3-position is the carbon adjacent to ring junction atoms in ring A
    beta_oh_pattern = Chem.MolFromSmarts('''
        [C@H]([O])[#6]
    ''')
    matches = mol.GetSubstructMatches(beta_oh_pattern)
    found_beta_oh = False
    for match in matches:
        atom_idx = match[0]  # The carbon with OH
        atom = mol.GetAtomWithIdx(atom_idx)
        # Check if the atom is in the steroid core
        if mol.GetSubstructMatch(steroid_pattern, useChirality=False):
            found_beta_oh = True
            break
    if not found_beta_oh:
        return False, "No 3beta-hydroxy group found"
    
    return True, "Contains steroid backbone with 3beta-hydroxy group"