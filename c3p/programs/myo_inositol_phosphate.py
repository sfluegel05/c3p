"""
Classifies: CHEBI:25448 myo-inositol phosphate
"""
"""
Classifies: myo-inositol phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_myo_inositol_phosphate(smiles: str):
    """
    Determines if a molecule is a myo-inositol phosphate based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a myo-inositol phosphate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count carbons (should be exactly 6 for myo-inositol)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 6:
        return False, f"Found {c_count} carbons, need exactly 6 for myo-inositol"
    
    # Check for phosphorus
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    if p_count == 0:
        return False, "No phosphorus atoms found"
    
    # Check for ring size (should be 6)
    ring_info = mol.GetRingInfo()
    ring_sizes = [len(ring) for ring in ring_info.AtomRings()]
    if not any(size == 6 for size in ring_sizes):
        return False, "No 6-membered ring found"
        
    # Multiple patterns to match different representations of myo-inositol core
    myo_patterns = [
        # Pattern 1: Standard myo-inositol with explicit stereochemistry
        "[C@@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)1",
        # Pattern 2: Alternative representation
        "[C@H]1[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)1",
        # Pattern 3: More general pattern with any substituents
        "C1C(O)C(O)C(O)C(O)C(O)1"
    ]
    
    found_myo = False
    for pattern in myo_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            found_myo = True
            break
            
    if not found_myo:
        return False, "Does not match myo-inositol core pattern"
    
    # Check for phosphate groups - multiple patterns to match different forms
    phosphate_patterns = [
        # Standard phosphate
        "[OX2][P]([OX2,OX1-])([OX2,OX1-])=[OX1,OX2-]",
        # Charged phosphate
        "[O-][P](=[O])([O-])[O-]",
        # Alternative phosphate representation
        "[P](=O)([O,O-])([O,O-])[O,O-]"
    ]
    
    phosphate_found = False
    for pattern in phosphate_patterns:
        phos_pattern = Chem.MolFromSmarts(pattern)
        if mol.HasSubstructMatch(phos_pattern):
            phosphate_found = True
            break
            
    if not phosphate_found:
        return False, "No phosphate groups found"
    
    # Verify phosphate attachment to ring
    ring_atoms = ring_info.AtomRings()[0]
    ring_carbons = set(idx for idx in ring_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
    
    phosphate_attached = False
    for carbon_idx in ring_carbons:
        carbon = mol.GetAtomWithIdx(carbon_idx)
        for neighbor in carbon.GetNeighbors():
            if neighbor.GetAtomicNum() == 8:  # Oxygen
                for next_neighbor in neighbor.GetNeighbors():
                    if next_neighbor.GetAtomicNum() == 15:  # Phosphorus
                        phosphate_attached = True
                        break
                if phosphate_attached:
                    break
        if phosphate_attached:
            break
            
    if not phosphate_attached:
        return False, "Phosphate groups not properly attached to inositol ring"

    # Count total oxygens to verify complete phosphate groups
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    min_o_needed = 6 + (p_count * 3)  # 6 for ring hydroxyls + 3 per phosphate
    if o_count < min_o_needed:
        return False, "Insufficient oxygen atoms for complete phosphate groups"

    return True, f"Contains myo-inositol core with {p_count} phosphorus atoms in phosphate groups"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:25448',
                          'name': 'myo-inositol phosphate',
                          'definition': 'An inositol phosphate in which the '
                                        'inositol component has '
                                        'myo-configuration.',
                          'parents': ['CHEBI:24846'],
                          'xrefs': ['PMID:1694860'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.6929133858267718 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: '
               'O[C@@H]1[C@@H](OP(O)(O)=O)[C@H](OP(O)(O)=O)[C@@H](OP(O)(O)=O)[C@H](OP(O)(O)=O)[C@@H]1OP(O)(O)=O '
               'NAME: 1D-myo-inositol 1,2,4,5,6-pentakisphosphate REASON: '
               'CORRECT Contains myo-inositol core with 5 phosphorus atoms in '
               'phosphate groups\n'
               ' * SMILES: '
               'OP(O)(=O)O[C@H]1[C@H](OP(O)(O)=O)[C@@H](OP(O)(=O)OP(O)(O)=O)[C@@H](OP(O)(O)=O)[C@@H](OP(O)(=O)OP(O)(O)=O)[C@@H]1OP(O)(O)=O '
               'NAME: 1,3-bis(diphospho)-myo-inositol '
               '2,4,5,6-tetrakisphosphate REASON: CORRECT Contains '
               'myo-inositol core with 8 phosphorus atoms in phosphate groups\n'
               ' * SMILES: '
               'OP(O)(=O)O[C@H]1[C@H](OP(O)(O)=O)[C@@H](OP(O)(O)=O)[C@H](OP(O)(=O)OP(O)(O)=O)[C@H](OP(O)(O)=O)[C@@H]1OP(O)(O)=O '
               'NAME: 1-diphospho-1D-myo-inositol 2,3,4,5,6-pentakisphosphate '
               'REASON: CORRECT Contains myo-inositol core with 7 phosphorus '
               'atoms in phosphate groups\n'
               ' * SMILES: '
               'O[C@@H]1[C@@H](O)[C@H](O)[C@@H](OP(O)(O)=O)[C@H](O)[C@@H]1O '
               'NAME: 1D-myo-inositol 6-phosphate REASON: CORRECT Contains '
               'myo-inositol core with 1 phosphorus atoms in phosphate groups\n'
               ' * SMILES: '
               'OP(O)(=O)O[C@@H]1[C@@H](OP(O)(O)=O)[C@H](OP(O)(=O)OP(O)(O)=O)[C@@H](OP(O)(O)=O)[C@H](OP(O)(=O)OP(O)(O)=O)[C@@H]1OP(O)(O)=O '
               'NAME: 1,5-bis(diphospho)-1D-myo-inositol '
               '2,3,4,6-tetrakisphosphate REASON: CORRECT Contains '
               'myo-inositol core with 8 phosphorus atoms in phosphate groups\n'
               ' * SMILES: '
               '[C@H]1([C@@H]([C@@H]([C@H]([C@@H]([C@H]1O)O)OP(O)(=O)O)O)OP(=O)(O)O)O '
               'NAME: 1D-myo-inositol 2,4-bisphosphate REASON: CORRECT '
               'Contains myo-inositol core with 2 phosphorus atoms in '
               'phosphate groups\n'
               ' * SMILES: '
               'O[C@H]1[C@@H](O)[C@H](O)[C@@H](OP(O)(O)=O)[C@H](OP(O)(O)=O)[C@H]1O '
               'NAME: 1D-myo-inositol 4,5-bisphosphate REASON: CORRECT '
               'Contains myo-inositol core with 2 phosphorus atoms in '
               'phosphate groups\n'
               ' * SMILES: '
               'O[C@H]1[C@H](O)[C@@H](OP(O)(O)=O)[C@@H](O)[C@@H](OP(O)(O)=O)[C@@H]1O '
               'NAME: myo-inositol 1,3-bisphosphate REASON: CORRECT Contains '
               'myo-inositol core with 2 phosphorus atoms in phosphate groups\n'
               ' * SMILES: '
               'OP(O)(=O)O[C@@H]1[C@H](OP(O)(O)=O)[C@@H](OP(O)(O)=O)[C@H](OP(O)(=O)OP(O)(=O)OP(O)(O)=O)[C@@H](OP(O)(O)=O)[C@@H]1OP(O)(O)=O '
               'NAME: 5-triphospho-1D-myo-inositol 1,2,3,4,6-pentakisphosphate '
               'REASON: CORRECT Contains myo-inositol core with 8 phosphorus '
               'atoms in phosphate groups\n'
               ' * SMILES: '
               'O[C@@H]1[C@H](OP(O)(O)=O)[C@@H](OP(O)(O)=O)[C@H](OP(O)(O)=O)[C@@H](OP(O)(O)=O)[C@@H]1OP(O)(=O)OP(O)(O)=O '
               'NAME: 1-diphospho-1D-myo-inositol 3,4,5,6-tetrakisphosphate '
               'REASON: CORRECT Contains myo-inositol core with 6 phosphorus '
               'atoms in phosphate groups\n'
               ' * SMILES: '
               '[C@@H]1(O)[C@H](OP([O-])(=O)[O-])[C@H]([C@H](OP(OP([O-])(=O)[O-])(=O)[O-])[C@H]([C@@H]1OP(=O)([O-])[O-])OP([O-])(=O)[O-])OP([O-])(=O)[O-] '
               'NAME: 5-diphospho-1D-myo-inositol '
               '1,3,4,6-tetrakisphosphate(11-) REASON: CORRECT Contains '
               'myo-inositol core with 6 phosphorus atoms in phosphate groups\n'
               ' * SMILES: '
               'O[C@H]1[C@H](OP(O)(O)=O)[C@@H](OP(O)(O)=O)[C@@H](OP(O)(O)=O)[C@@H](OP(O)(O)=O)[C@@H]1OP(O)(O)=O '
               'NAME: 1L-myo-inositol 1,2,3,4,6-pentakisphosphate REASON: '
               'CORRECT Contains myo-inositol core with 5 phosphorus atoms in '
               'phosphate groups\n'
               ' * SMILES: '
               '[C@@H]1(OP(=O)(O)O)[C@@H](OP(=O)(O)OP(O)(O)=O)[C@H](OP(=O)(O)O)[C@@H](OP(=O)(O)OP(O)(O)=O)[C@@H](OP(=O)(O)O)[C@H]1OP(=O)(O)O '
               'NAME: 3,5-bis(diphospho)-1D-myo-inositol '
               '1,2,4,6-tetrakisphosphate REASON: CORRECT Contains '
               'myo-inositol core with 8 phosphorus atoms in phosphate groups\n'
               ' * SMILES: '
               'O[C@@H]1[C@H](O)[C@@H](OP(O)(O)=O)[C@H](OP(O)(O)=O)[C@@H](O)[C@@H]1OP(O)(O)=O '
               'NAME: 1D-myo-inositol 1,4,5-trisphosphate REASON: CORRECT '
               'Contains myo-inositol core with 3 phosphorus atoms in '
               'phosphate groups\n'
               ' * SMILES: '
               'O[C@H]1[C@H](O)[C@@H](O)[C@@H](OP(O)(O)=O)[C@@H](O)[C@@H]1O '
               'NAME: 1D-myo-inositol 2-phosphate REASON: CORRECT Contains '
               'myo-inositol core with 1 phosphorus atoms in phosphate groups\n'
               ' * SMILES: '
               'O[C@H]1[C@H](O)[C@@H](O)[C@H](OP(O)(O)=O)[C@@H](OP(O)(O)=O)[C@H]1O '
               'NAME: 1D-myo-inositol 3,4-bisphosphate REASON: CORRECT '
               'Contains myo-inositol core with 2 phosphorus atoms in '
               'phosphate groups\n'
               ' * SMILES: '
               'O[C@H]1[C@@H](OP(O)(O)=O)[C@H](OP(O)(O)=O)[C@@H](OP(O)(O)=O)[C@H](OP(O)(O)=O)[C@H]1OP(O)(=O)OP(O)(O)=O '
               'NAME: 3-diphospho-1D-myo-inositol 1,4,5,6-tetrakisphosphate '
               'REASON: CORRECT Contains myo-inositol core with 6 phosphorus '
               'atoms in phosphate groups\n'
               ' * SMILES: '
               'O[C@H]1[C@@H](OP(O)(O)=O)[C@H](O)[C@@H](OP(O)(O)=O)[C@H](OP(O)(O)=O)[C@H]1OP(O)(O)=O '
               'NAME: 1D-myo-inositol 1,3,4,5-tetrakisphosphate REASON: '
               'CORRECT Contains myo-inositol core with 4 phosphorus atoms in '
               'phosphate groups\n'
               ' * SMILES: '
               'OP(O)(=O)O[C@H]1[C@@H](OP(O)(O)=O)[C@H](OP(O)(O)=O)[C@@H](OP(O)(=O)OP(O)(O)=O)[C@H](OP(O)(=O)OP(O)(O)=O)[C@H]1OP(O)(O)=O '
               'NAME: 4,5-bis(diphospho)-1D-myo-inositol tetrakisphosphate '
               'REASON: CORRECT Contains myo-inositol core with 8 phosphorus '
               'atoms in phosphate groups\n'
               ' * SMILES: '
               'O[C@H]1[C@H](O)[C@H](O)[C@H](OP(O)(O)=O)[C@@H](O)[C@@H]1O '
               'NAME: 1D-myo-inositol 3-phosphate REASON: CORRECT Contains '
               'myo-inositol core with 1 phosphorus atoms in phosphate groups\n'
               ' * SMILES: '
               'OP(O)(=O)O[C@@H]1[C@H](OP(O)(O)=O)[C@@H](OP(O)(O)=O)[C@H](OP(O)(=O)OP(O)(O)=O)[C@@H](OP(O)(O)=O)[C@@H]1OP(O)(O)=O '
               'NAME: 5-diphospho-1D-myo-inositol pentakisphosphate REASON: '
               'CORRECT Contains myo-inositol core with 7 phosphorus atoms in '
               'phosphate groups\n'
               ' * SMILES: '
               'O[C@H]1[C@H](OP(O)(O)=O)[C@H](OP(O)(O)=O)[C@H](OP(O)(O)=O)[C@@H](OP(O)(O)=O)[C@@H]1OP(O)(O)=O '
               'NAME: 1D-myo-inositol 1,2,3,4,5-pentakisphosphate REASON: '
               'CORRECT Contains myo-inositol core with 5 phosphorus atoms in '
               'phosphate groups\n'
               ' * SMILES: '
               'O[C@H]1[C@H](OP(O)(O)=O)[C@H](O)[C@H](OP(O)(O)=O)[C@@H](O)[C@@H]1OP(O)(O)=O '
               'NAME: myo-inositol 1,3,5-trisphosphate REASON: CORRECT '
               'Contains myo-inositol core with 3 phosphorus atoms in '
               'phosphate groups\n'
               ' * SMILES: '
               'O[C@@H]1[C@H](O)[C@@H](O)[C@H](OP(O)(O)=O)[C@@H](O)[C@@H]1O '
               'NAME: 1D-myo-inositol 5-phosphate REASON: CORRECT Contains '
               'myo-inositol core with 1 phosphorus atoms in phosphate groups\n'
               ' * SMILES: '
               'O[C@H]1[C@@H](O)[C@H](OP(O)(O)=O)[C@@H](OP(O)(O)=O)[C@H](OP(O)(O)=O)[C@H]1OP(O)(O)=O '
               'NAME: 1D-myo-inositol 3,4,5,6-tetrakisphosphate REASON: '
               'CORRECT Contains myo-inositol core with 4 phosphorus atoms in '
               'phosphate groups\n'
               'False positives: SMILES: '
               '[C@@H]1(OP(=O)([O-])[O-])[C@@H](O)[C@@H]([C@H]([C@H]([C@H]1OP([O-])(=O)[O-])OP(=O)([O-])[O-])O)O '
               'NAME: 1D-myo-inositol 1,2,6-trisphosphate(6-) REASON: WRONGLY '
               'CLASSIFIED Contains myo-inositol core with 3 phosphorus atoms '
               'in phosphate groups\n'
               ' * SMILES: '
               'O[C@H]1[C@@H](OP([O-])([O-])=O)[C@H](O)[C@@H](OP([O-])([O-])=O)[C@H](OP([O-])([O-])=O)[C@H]1OP([O-])([O-])=O '
               'NAME: 1D-myo-inositol 1,3,4,5-tetrakisphosphate(8-) REASON: '
               'WRONGLY CLASSIFIED Contains myo-inositol core with 4 '
               'phosphorus atoms in phosphate groups\n'
               ' * SMILES: P(OC1[C@H](O)[C@@H](O)C(O)[C@H](O)[C@@H]1O)(O)(O)=O '
               'NAME: Inositol 4-phosphate REASON: WRONGLY CLASSIFIED Contains '
               'myo-inositol core with 1 phosphorus atoms in phosphate groups\n'
               ' * SMILES: '
               'P(OC1C(OP(O)(O)=O)C(OP(O)(O)=O)C(OP(O)(O)=O)C(OP(O)(O)=O)C1OP(O)(O)=O)(O)(O)=O '
               'NAME: Diphosphoinositol tetrakisphosphate REASON: WRONGLY '
               'CLASSIFIED Contains myo-inositol core with 6 phosphorus atoms '
               'in phosphate groups\n'
               ' * SMILES: '
               '[C@H]1([C@@H]([C@@H]([C@H]([C@@H]([C@H]1O)OP(=O)([O-])[O-])O)O)OP(=O)([O-])[O-])O '
               'NAME: 1D-myo-inositol 2,5-bisphosphate(4-) REASON: WRONGLY '
               'CLASSIFIED Contains myo-inositol core with 2 phosphorus atoms '
               'in phosphate groups\n'
               ' * SMILES: '
               'O[C@H]1[C@H](O)[C@@H](O)[C@H](OP([O-])([O-])=O)[C@H](O)[C@@H]1O '
               'NAME: 1D-myo-inositol 1-phosphate(2-) REASON: WRONGLY '
               'CLASSIFIED Contains myo-inositol core with 1 phosphorus atoms '
               'in phosphate groups\n'
               ' * SMILES: '
               '[C@@H]1(O)[C@@H](O)[C@@H]([C@H]([C@H]([C@H]1OP([O-])(=O)[O-])OP(=O)([O-])[O-])O)O '
               'NAME: 1D-myo-inositol 1,2-bisphosphate(4-) REASON: WRONGLY '
               'CLASSIFIED Contains myo-inositol core with 2 phosphorus atoms '
               'in phosphate groups\n'
               ' * SMILES: '
               'O[C@@H]1[C@H](O)[C@@H](OP([O-])([O-])=O)[C@H](OP([O-])([O-])=O)[C@@H](O)[C@@H]1OP([O-])([O-])=O '
               'NAME: 1D-myo-inositol 1,4,5-trisphosphate(6-) REASON: WRONGLY '
               'CLASSIFIED Contains myo-inositol core with 3 phosphorus atoms '
               'in phosphate groups\n'
               ' * SMILES: '
               '[C@H]1(OP([O-])(OP([O-])(=O)[O-])=O)[C@@H](OP([O-])(=O)[O-])[C@@H](OP([O-])(=O)[O-])[C@@H]([C@H]([C@@H]1OP([O-])(=O)[O-])OP(=O)([O-])[O-])O '
               'NAME: 3-diphospho-1D-myo-inositol '
               '1,2,4,5-tetrakisphosphate(11-) REASON: WRONGLY CLASSIFIED '
               'Contains myo-inositol core with 6 phosphorus atoms in '
               'phosphate groups\n'
               ' * SMILES: '
               'O[C@H]1[C@H](OP([O-])([O-])=O)[C@@H](OP([O-])([O-])=O)[C@H](OP([O-])([O-])=O)[C@H](OP([O-])([O-])=O)[C@@H]1OP([O-])([O-])=O '
               'NAME: 1D-myo-inositol 1,2,3,5,6-pentakisphosphate(10-) REASON: '
               'WRONGLY CLASSIFIED Contains myo-inositol core with 5 '
               'phosphorus atoms in phosphate groups\n'
               ' * SMILES: '
               '[C@H]1([C@@H]([C@@H]([C@H]([C@@H]([C@H]1O)OP(=O)([O-])[O-])OP(=O)([O-])[O-])O)OP(=O)([O-])[O-])O '
               'NAME: 1D-myo-inositol 2,4,5-triphosphate(6-) REASON: WRONGLY '
               'CLASSIFIED Contains myo-inositol core with 3 phosphorus atoms '
               'in phosphate groups\n'
               ' * SMILES: '
               '[C@H]1([C@@H]([C@@H]([C@H]([C@@H]([C@H]1O)OP(=O)([O-])[O-])OP(=O)([O-])[O-])OP([O-])(=O)[O-])OP(=O)([O-])[O-])O '
               'NAME: 1D-myo-inositol 2,3,4,5-tetrakisphosphate(8-) REASON: '
               'WRONGLY CLASSIFIED Contains myo-inositol core with 4 '
               'phosphorus atoms in phosphate groups\n'
               ' * SMILES: '
               '[O-]P([O-])(=O)O[C@@H]1[C@@H](OP([O-])([O-])=O)[C@H](OP([O-])(=O)OP([O-])([O-])=O)[C@@H](OP([O-])([O-])=O)[C@H](OP([O-])(=O)OP([O-])([O-])=O)[C@@H]1OP([O-])([O-])=O '
               'NAME: 1,5-bis(diphospho)-1D-myo-inositol '
               '2,3,4,6-tetrakisphosphate(14-) REASON: WRONGLY CLASSIFIED '
               'Contains myo-inositol core with 8 phosphorus atoms in '
               'phosphate groups\n'
               ' * SMILES: '
               '[C@@H]1(OP(=O)([O-])[O-])[C@@H](OP(=O)([O-])[O-])[C@H](OP(=O)([O-])[O-])[C@@H](OP(=O)([O-])OP([O-])([O-])=O)[C@@H](OP(=O)([O-])[O-])[C@H]1OP(=O)([O-])[O-] '
               'NAME: 3-diphospho-1D-myo-inositol '
               '1,2,4,5,6-pentakisphosphate(13-) REASON: WRONGLY CLASSIFIED '
               'Contains myo-inositol core with 7 phosphorus atoms in '
               'phosphate groups\n'
               ' * SMILES: '
               '[C@@H]1([C@@H]([C@@H]([C@@H]([C@H]([C@@H]1O)OP([O-])(=O)[O-])OP(=O)([O-])[O-])O)O)OP(=O)([O-])[O-] '
               'NAME: 1D-myo-inositol 3,4,6-trisphosphate(6-) REASON: WRONGLY '
               'CLASSIFIED Contains myo-inositol core with 3 phosphorus atoms '
               'in phosphate groups\n'
               ' * SMILES: '
               'O[C@H]1[C@H](O)[C@@H](O)[C@H]2OP(O)(=O)O[C@H]2[C@@H]1O NAME: '
               '1D-myo-inositol 1,2-cyclic phosphate REASON: WRONGLY '
               'CLASSIFIED Contains myo-inositol core with 1 phosphorus atoms '
               'in phosphate groups\n'
               ' * SMILES: '
               '[C@@H]1(OP(=O)([O-])[O-])[C@@H](OP(=O)([O-])OP([O-])([O-])=O)[C@H](OP(=O)([O-])[O-])[C@@H](OP(=O)([O-])OP([O-])(O)=O)[C@@H](OP(=O)([O-])[O-])[C@H]1OP(=O)([O-])[O-] '
               'NAME: 3,5-bis(diphospho)-1D-myo-inositol '
               '1,2,4,6-tetrakisphosphate(13-) REASON: WRONGLY CLASSIFIED '
               'Contains myo-inositol core with 8 phosphorus atoms in '
               'phosphate groups\n'
               ' * SMILES: '
               'O[C@H]1[C@@H](O)[C@@H](OP([O-])([O-])=O)[C@H](OP([O-])([O-])=O)[C@@H](OP([O-])([O-])=O)[C@@H]1OP([O-])([O-])=O '
               'NAME: 1D-myo-inositol 1,4,5,6-tetrakisphosphate(8-) REASON: '
               'WRONGLY CLASSIFIED Contains myo-inositol core with 4 '
               'phosphorus atoms in phosphate groups\n'
               ' * SMILES: '
               '[O-]P([O-])(=O)O[C@@H]1[C@H](OP([O-])([O-])=O)[C@@H](OP([O-])([O-])=O)[C@H](OP([O-])(=O)OP([O-])([O-])=O)[C@@H](OP([O-])([O-])=O)[C@@H]1OP([O-])([O-])=O '
               'NAME: 5-diphospho-1D-myo-inositol pentakisphosphate(13-) '
               'REASON: WRONGLY CLASSIFIED Contains myo-inositol core with 7 '
               'phosphorus atoms in phosphate groups\n'
               ' * SMILES: '
               'O[C@@H]1[C@H](O)[C@H](OP([O-])([O-])=O)[C@@H](O)[C@H](OP([O-])([O-])=O)[C@H]1O '
               'NAME: 1D-myo-inositol 3,5-bisphosphate(4-) REASON: WRONGLY '
               'CLASSIFIED Contains myo-inositol core with 2 phosphorus atoms '
               'in phosphate groups\n'
               ' * SMILES: '
               'O[C@@H]1[C@H](O)[C@H](O)[C@@H](OP([O-])([O-])=O)[C@H](O)[C@H]1O '
               'NAME: 1D-myo-inositol 4-phosphate(2-) REASON: WRONGLY '
               'CLASSIFIED Contains myo-inositol core with 1 phosphorus atoms '
               'in phosphate groups\n'
               ' * SMILES: '
               'O[C@@H]1[C@@H](OP([O-])([O-])=O)[C@H](OP([O-])([O-])=O)[C@@H](OP([O-])([O-])=O)[C@H](OP([O-])([O-])=O)[C@@H]1OP([O-])([O-])=O '
               'NAME: 1D-myo-inositol 1,2,4,5,6-pentakisphosphate(10-) REASON: '
               'WRONGLY CLASSIFIED Contains myo-inositol core with 5 '
               'phosphorus atoms in phosphate groups\n'
               ' * SMILES: '
               '[C@@H]1([C@@H]([C@H]([C@@H]([C@H]([C@@H]1O)OP(=O)([O-])[O-])O)O)O)O '
               'NAME: scyllo-inositol 1-phosphate(2-) REASON: WRONGLY '
               'CLASSIFIED Contains myo-inositol core with 1 phosphorus atoms '
               'in phosphate groups\n'
               ' * SMILES: '
               'O[C@H]1[C@H](O)[C@@H](O)[C@H]2OP([O-])(=O)O[C@H]2[C@@H]1O '
               'NAME: 1D-myo-inositol 1,2-cyclic phosphate(1-) REASON: WRONGLY '
               'CLASSIFIED Contains myo-inositol core with 1 phosphorus atoms '
               'in phosphate groups\n'
               ' * SMILES: '
               'O[C@H]1[C@H](OP([O-])([O-])=O)[C@@H](OP([O-])([O-])=O)[C@@H](O)[C@@H](OP([O-])([O-])=O)[C@@H]1OP([O-])([O-])=O '
               'NAME: myo-inositol 1,3,4,6-tetrakisphosphate(8-) REASON: '
               'WRONGLY CLASSIFIED Contains myo-inositol core with 4 '
               'phosphorus atoms in phosphate groups\n'
               'False negatives: NONE\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'C[C@H]1CN(C(=O)C2=C(C3=CC=CC=C3CO[C@@H]1CN(C)C(=O)NC4=CC=C(C=C4)OC)C5=CC=CC=C5N2C)[C@H](C)CO',
                                     'name': 'LSM-30659',
                                     'reason': 'Found 34 carbons, need exactly '
                                               '6 for myo-inositol'},
                                 {   'smiles': 'c1c[nH+]c[nH]1',
                                     'name': 'imidazolium cation',
                                     'reason': 'Found 3 carbons, need exactly '
                                               '6 for myo-inositol'},
                                 {   'smiles': 'ClC=1C(O)=C2O[C@@]([C@H](NC(=O)[C@@H](NC(=O)[C@@H](NC)[C@@H](O)C(=C2)C1)C(C)=C)C(=O)N3[C@@H](C=CC3)C(=O)N/C(=C(/CC)\\C)/C(=O)N/C(=C/C(O)=O)/C(O)=O)(CC)C',
                                     'name': 'Phomopsin A',
                                     'reason': 'Found 36 carbons, need exactly '
                                               '6 for myo-inositol'},
                                 {   'smiles': 'ClC(Cl)[C@H](O)CC=1OC(=O)C=2C(O)=CC(=CC2C1)O',
                                     'name': 'Desmethyldichlorodiaportin',
                                     'reason': 'Found 12 carbons, need exactly '
                                               '6 for myo-inositol'},
                                 {   'smiles': 'ClC=1C(=C(O)C2=C(C1O)C(=O)C=CC2=O)CC=C(C)C',
                                     'name': 'Chlorosesamone',
                                     'reason': 'Found 15 carbons, need exactly '
                                               '6 for myo-inositol'},
                                 {   'smiles': 'CCCCC[C@H]1O[C@H]1C\\C=C/CCCCCCCC(O)=O',
                                     'name': '(+)-vernolic acid',
                                     'reason': 'Found 18 carbons, need exactly '
                                               '6 for myo-inositol'},
                                 {   'smiles': 'O(C1C(O)C(OC(OC=2C3=C(C=CC2)C=C(C(=C3O)C(=O)C)C)C1O)CO)C4OC(C(O)C(O)C4O)CO',
                                     'name': 'Orientaloside',
                                     'reason': 'Found 25 carbons, need exactly '
                                               '6 for myo-inositol'},
                                 {   'smiles': 'O(C1=C(C=CC=C1C(O)=O)C)C(=O)C',
                                     'name': 'CRESOPYRINE',
                                     'reason': 'Found 10 carbons, need exactly '
                                               '6 for myo-inositol'},
                                 {   'smiles': 'C[C@H]1CN(S(=O)(=O)C2=C(C=C(C=C2)C#CC3=CC=NC=C3)O[C@@H]1CN(C)C(=O)NC4=CC5=C(C=C4)OCO5)[C@@H](C)CO',
                                     'name': '3-(1,3-benzodioxol-5-yl)-1-[[(4S,5S)-2-[(2S)-1-hydroxypropan-2-yl]-4-methyl-1,1-dioxo-8-(2-pyridin-4-ylethynyl)-4,5-dihydro-3H-6,1$l^{6},2-benzoxathiazocin-5-yl]methyl]-1-methylurea',
                                     'reason': 'Found 30 carbons, need exactly '
                                               '6 for myo-inositol'},
                                 {   'smiles': 'O1C2C3C(CCC3=C)C(CCC2C(C1=O)=C)=C',
                                     'name': '3,6,9-Trimethylidene-3a,4,5,6a,7,8,9a,9b-octahydroazuleno[4,5-b]furan-2-one',
                                     'reason': 'Found 15 carbons, need exactly '
                                               '6 for myo-inositol'}],
    'sample_false_negatives': [],
    'attempt': 3,
    'reasoning': 'Looking at the outcomes, I notice that the main issue is '
                 'that the program is incorrectly handling charged forms of '
                 'myo-inositol phosphates. The false positives are mostly '
                 'charged versions of valid myo-inositol phosphates (with [O-] '
                 'groups), but the program should handle both charged and '
                 'uncharged forms equally.\n'
                 '\n'
                 "Here's my analysis of what needs to be improved:\n"
                 '\n'
                 "1. The uncharge_smiles step isn't effectively normalizing "
                 'the structures - we should either properly handle both '
                 'charged and uncharged forms, or normalize them more '
                 'effectively.\n'
                 '\n'
                 '2. The myo-inositol stereochemistry pattern might be too '
                 'strict or not correctly matching all valid configurations.\n'
                 '\n'
                 '3. We need to better handle the various phosphate group '
                 'patterns that can appear in both charged and uncharged '
                 'forms.\n'
                 '\n'
                 "Here's my improved version:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 44,
    'num_false_positives': 36,
    'num_true_negatives': 142220,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.55,
    'recall': 1.0,
    'f1': 0.7096774193548387,
    'accuracy': 0.999747013352073,
    'negative_predictive_value': 1.0}