"""
Classifies: CHEBI:140331 4'-hydroxyflavanones
"""
"""
Classifies: 4'-hydroxyflavanones
Definition: Any hydroxyflavanone having a hydroxy substituent located at position 4'.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_4__hydroxyflavanones(smiles: str):
    """
    Determines if a molecule is a 4'-hydroxyflavanone based on its SMILES string.
    A 4'-hydroxyflavanone has a flavanone backbone with a hydroxy group at the 4' position of the phenyl ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 4'-hydroxyflavanone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more flexible flavanone backbone pattern
    flavanone_pattern = Chem.MolFromSmarts("[#6]1-[#6](-[#6]=O)-[#6](-[#8]-[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2)-[#6]-1")
    if not mol.HasSubstructMatch(flavanone_pattern):
        return False, "No flavanone backbone found"

    # Find the phenyl ring attached to the flavanone core
    phenyl_pattern = Chem.MolFromSmarts("[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1")
    phenyl_matches = mol.GetSubstructMatches(phenyl_pattern)
    
    # Check each phenyl ring for a hydroxy group at the 4' position
    for match in phenyl_matches:
        # Get the atoms in the phenyl ring
        ring_atoms = [mol.GetAtomWithIdx(i) for i in match]
        
        # Check if this phenyl ring is attached to the flavanone oxygen
        for atom in ring_atoms:
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'O' and neighbor.GetDegree() == 2:
                    # This is the phenyl ring attached to the flavanone core
                    # Now check for a hydroxy group at the 4' position
                    # The 4' position is the carbon opposite the attachment point
                    attachment_point = atom.GetIdx()
                    opposite_carbon = None
                    for ring_atom in ring_atoms:
                        if ring_atom.GetIdx() != attachment_point and \
                           not mol.GetBondBetweenAtoms(ring_atom.GetIdx(), attachment_point):
                            opposite_carbon = ring_atom
                            break
                    
                    if opposite_carbon:
                        # Check if the opposite carbon has a hydroxy group
                        for neighbor in opposite_carbon.GetNeighbors():
                            if neighbor.GetSymbol() == 'O' and neighbor.GetTotalNumHs() == 1:
                                return True, "Contains flavanone backbone with a hydroxy group at the 4' position"
                        # Also check for deprotonated hydroxy groups
                        for neighbor in opposite_carbon.GetNeighbors():
                            if neighbor.GetSymbol() == 'O' and neighbor.GetFormalCharge() == -1:
                                return True, "Contains flavanone backbone with a deprotonated hydroxy group at the 4' position"

    return False, "No hydroxy group at the 4' position"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:XXXXX',
                          'name': '4\'-hydroxyflavanones',
                          'definition': 'Any hydroxyflavanone having a hydroxy substituent located at position 4\'.',
                          'parents': ['CHEBI:XXXXX', 'CHEBI:XXXXX']},
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
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': 0.0}