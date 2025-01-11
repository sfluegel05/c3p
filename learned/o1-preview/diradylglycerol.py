"""
Classifies: CHEBI:76578 diradylglycerol
"""
"""
Classifies: CHEBI:18185 diradylglycerol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_diradylglycerol(smiles: str):
    """
    Determines if a molecule is a diradylglycerol based on its SMILES string.
    A diradylglycerol is a glycerol molecule with two of its three hydroxyl groups
    substituted by acyl, alkyl, or alk-1-enyl groups at any positions.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a diradylglycerol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string to RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define glycerol backbone pattern
    # Glycerol: C-C-C with hydroxyl groups attached to each carbon
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)CO")
    matches = mol.GetSubstructMatches(glycerol_pattern)
    if not matches:
        return False, "No glycerol backbone found"
    
    # Find glycerol carbon atoms
    glycerol_atoms = set()
    for match in matches:
        glycerol_atoms.update(match)
    
    # Check substitutions on glycerol carbons
    substitution_count = 0
    for atom_idx in glycerol_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        # Count non-hydrogen neighbors excluding the glycerol oxygen atoms
        neighbor_atoms = [nbr for nbr in atom.GetNeighbors() if nbr.GetIdx() not in glycerol_atoms]
        for nbr in neighbor_atoms:
            # Check if neighbor is connected via ester (acyl), ether (alkyl), or vinyl ether (alk-1-enyl)
            bond = mol.GetBondBetweenAtoms(atom_idx, nbr.GetIdx())
            bond_type = bond.GetBondType()
            nbr_atom_num = nbr.GetAtomicNum()
            if bond_type == Chem.BondType.SINGLE:
                if nbr_atom_num == 6:
                    # Possible alkyl or alk-1-enyl group
                    # Check for further unsaturation or ester linkage
                    # Acyl group would be identified by a carbonyl group adjacent to the oxygen
                    # For simplicity, allow any carbon chain (alkyl)
                    substitution_count += 1
                elif nbr_atom_num == 8:
                    # Possible ester linkage (acyl group)
                    # Check for carbonyl adjacent to the oxygen
                    ester_pattern = Chem.MolFromSmarts("[$(O=C)]([#6])[O]")
                    if mol.HasSubstructMatch(ester_pattern):
                        substitution_count += 1
            elif bond_type == Chem.BondType.DOUBLE:
                # Possible alk-1-enyl group
                substitution_count += 1
    
    if substitution_count != 2:
        return False, f"Found {substitution_count} substitutions, need exactly 2"
    
    return True, "Contains glycerol backbone with two substituents (acyl, alkyl, or alk-1-enyl)"

__metadata__ = {   'chemical_class': {   'id': 'CHEBI:18185',
                              'name': 'diradylglycerol',
                              'definition': 'Any lipid that is glycerol bearing '
                                            'two substituent groups - either acyl, '
                                            'alkyl, or alk-1-enyl - at any two of '
                                            'the three possible positions.',
                              'parents': []},
    'config': {   'llm_model_name': 'lbl/claude-sonnet',
                  'f1_threshold': 0.8,
                  'max_attempts': 5,
                  'max_positive_instances': None,
                  'max_positive_to_test': None,
                  'max_negative_to_test': None,
                  'max_positive_in_prompt': 50,
                  'max_negative_in_prompt': 20,
                  'max_instances_in_prompt': 100,
                  'test_proportion': 0.1},
    'message': None,
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': None,
    'num_false_positives': None,
    'num_true_negatives': None,
    'num_false_negatives': None,
    'num_negatives': None,
    'precision': None,
    'recall': None,
    'f1': None,
    'accuracy': None}