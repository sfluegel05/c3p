"""
Classifies: CHEBI:76579 triradylglycerol
"""
"""
Classifies: CHEBI:36085 triradylglycerol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_triradylglycerol(smiles: str):
    """
    Determines if a molecule is a triradylglycerol based on its SMILES string.
    A triradylglycerol is a glycerol backbone with three substituent groups (acyl, alkyl, or alk-1-enyl)
    attached via oxygen atoms at positions sn-1, sn-2, and sn-3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triradylglycerol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define substituted glycerol backbone pattern
    glycerol_pattern = Chem.MolFromSmarts("[C:1]-[C:2]-[C:3]")
    
    # Find glycerol backbone matches
    backbone_matches = mol.GetSubstructMatches(glycerol_pattern)
    if not backbone_matches:
        return False, "No glycerol backbone found"

    for backbone_match in backbone_matches:
        c1_idx, c2_idx, c3_idx = backbone_match

        # Check that each carbon is connected to an oxygen
        c1 = mol.GetAtomWithIdx(c1_idx)
        c2 = mol.GetAtomWithIdx(c2_idx)
        c3 = mol.GetAtomWithIdx(c3_idx)

        # Get oxygens attached to each carbon
        o1 = get_attached_oxygen(c1, exclude_idx=c2_idx)
        o2 = get_attached_oxygen(c2, exclude_idxs=[c1_idx, c3_idx])
        o3 = get_attached_oxygen(c3, exclude_idx=c2_idx)

        if None in (o1, o2, o3):
            continue  # One of the carbons doesn't have an attached oxygen

        # Now, check the substituents attached to oxygens
        substituents = []
        for i, (o_atom, c_idx) in enumerate([(o1, c1_idx), (o2, c2_idx), (o3, c3_idx)], 1):
            # Get the atom connected to oxygen (substituent)
            substituent = get_substituent(o_atom, exclude_idx=c_idx)
            if substituent is None:
                return False, f"No substituent attached to oxygen at position {i}"

            # Analyze the substituent
            if is_acyl_group(o_atom, substituent):
                substituents.append('acyl')
            elif is_alkyl_group(o_atom, substituent):
                substituents.append('alkyl')
            elif is_alk1enyl_group(o_atom, substituent):
                substituents.append('alk-1-enyl')
            else:
                return False, f"Substituent at position {i} is not acyl, alkyl, or alk-1-enyl group"

        # If all substituents are valid, return True
        return True, f"Contains glycerol backbone with substituents: {', '.join(substituents)}"

    return False, "No suitable glycerol backbone with the required substituents found"

def get_attached_oxygen(carbon_atom, exclude_idx=None, exclude_idxs=None):
    """
    Returns the oxygen atom attached to the given carbon atom.

    Args:
        carbon_atom: RDKit Atom object for carbon
        exclude_idx: Atom index to exclude from consideration
        exclude_idxs: List of atom indices to exclude

    Returns:
        RDKit Atom object for oxygen, or None if not found
    """
    if exclude_idxs is None:
        exclude_idxs = []
    if exclude_idx is not None:
        exclude_idxs.append(exclude_idx)

    for neighbor in carbon_atom.GetNeighbors():
        if neighbor.GetAtomicNum() == 8 and neighbor.GetIdx() not in exclude_idxs:
            return neighbor
    return None

def get_substituent(o_atom, exclude_idx):
    """
    Returns the substituent atom connected to the oxygen atom.

    Args:
        o_atom: RDKit Atom object for oxygen
        exclude_idx: Atom index to exclude (the carbon of the glycerol backbone)

    Returns:
        RDKit Atom object for substituent, or None if not found
    """
    for neighbor in o_atom.GetNeighbors():
        if neighbor.GetIdx() != exclude_idx:
            return neighbor
    return None

def is_acyl_group(o_atom, substituent):
    """
    Checks if the substituent attached to the oxygen atom is an acyl group.

    Args:
        o_atom: Oxygen atom (RDKit Atom object)
        substituent: Substituent atom connected to oxygen (RDKit Atom object)

    Returns:
        bool: True if acyl group, False otherwise
    """
    if substituent.GetAtomicNum() != 6:
        return False

    # Check for carbonyl group (C=O) attached to substituent
    has_carbonyl = False
    for bond in substituent.GetBonds():
        bond_type = bond.GetBondType()
        other_atom = bond.GetOtherAtom(substituent)
        if other_atom.GetIdx() == o_atom.GetIdx():
            continue  # Skip bond back to oxygen
        if bond_type == Chem.rdchem.BondType.DOUBLE and other_atom.GetAtomicNum() == 8:
            has_carbonyl = True
            break
    return has_carbonyl

def is_alkyl_group(o_atom, substituent):
    """
    Checks if the substituent attached to the oxygen atom is an alkyl group.

    Args:
        o_atom: Oxygen atom (RDKit Atom object)
        substituent: Substituent atom connected to oxygen (RDKit Atom object)

    Returns:
        bool: True if alkyl group, False otherwise
    """
    if substituent.GetAtomicNum() != 6:
        return False

    # Check that substituent is sp3 hybridized carbon with only single bonds
    if substituent.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
        return False
    for bond in substituent.GetBonds():
        if bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
            return False
    return True

def is_alk1enyl_group(o_atom, substituent):
    """
    Checks if the substituent attached to the oxygen atom is an alk-1-enyl group.

    Args:
        o_atom: Oxygen atom (RDKit Atom object)
        substituent: Substituent atom connected to oxygen (RDKit Atom object)

    Returns:
        bool: True if alk-1-enyl group, False otherwise
    """
    if substituent.GetAtomicNum() != 6:
        return False

    # Check for double bond between substituent and next carbon
    has_double_bond = False
    for bond in substituent.GetBonds():
        other_atom = bond.GetOtherAtom(substituent)
        if other_atom.GetIdx() == o_atom.GetIdx():
            continue  # Skip bond back to oxygen
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and other_atom.GetAtomicNum() == 6:
            has_double_bond = True
            break
    return has_double_bond

__metadata__ = {   'chemical_class': {   'id': 'CHEBI:36085',
                          'name': 'triradylglycerol',
                          'definition': 'A glycerol compound having one of three possible substituent groups - either acyl, alkyl, or alk-1-enyl - at each of the three possible positions sn-1, sn-2 or sn-3.',
                          'parents': ['CHEBI:35741', 'CHEBI:17754']},
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
    'attempt': 1,
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