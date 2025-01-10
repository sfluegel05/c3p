"""
Classifies: CHEBI:16337 phosphatidic acid
"""
"""
Classifies: phosphatidic acid
"""
from rdkit import Chem

def is_phosphatidic_acid(smiles: str):
    """
    Determines if a molecule is a phosphatidic acid based on its SMILES string.
    A phosphatidic acid is a derivative of glycerol in which one hydroxy group is esterified with phosphoric acid
    and the other two are esterified with fatty acids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidic acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone (C-C-C with 3 oxygens attached)
    # Using a more general pattern to account for stereochemistry and variances
    glycerol_pattern = Chem.MolFromSmarts("C([O])[C@H]([O])[C]([O])")
    matches = mol.GetSubstructMatches(glycerol_pattern)
    if not matches:
        # Try without stereochemistry
        glycerol_pattern = Chem.MolFromSmarts("C([O])C([O])C([O])")
        matches = mol.GetSubstructMatches(glycerol_pattern)
        if not matches:
            return False, "No glycerol backbone found"

    # For each glycerol backbone found, check attachments
    for match in matches:
        c1_idx, c2_idx, c3_idx, o1_idx, o2_idx, o3_idx = match[0], match[1], match[2], match[3], match[4], match[5]
        c1 = mol.GetAtomWithIdx(c1_idx)
        c2 = mol.GetAtomWithIdx(c2_idx)
        c3 = mol.GetAtomWithIdx(c3_idx)
        o1 = mol.GetAtomWithIdx(o1_idx)
        o2 = mol.GetAtomWithIdx(o2_idx)
        o3 = mol.GetAtomWithIdx(o3_idx)
        # Check for phosphate group attached to one oxygen
        phosphate_attached = False
        ester_count = 0
        for oxygen in [o1, o2, o3]:
            for neighbor in oxygen.GetNeighbors():
                if neighbor.GetAtomicNum() == 15:  # Phosphorus
                    phosphate_attached = True
                    break
                elif neighbor.GetAtomicNum() == 6:
                    # Check if connected to carbonyl (ester linkage)
                    for n_bond in neighbor.GetBonds():
                        if n_bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and n_bond.GetOtherAtom(neighbor).GetAtomicNum() == 8:
                            ester_count += 1
                            break
            if phosphate_attached:
                break  # Only one phosphate group needed

        if not phosphate_attached:
            return False, "No phosphate group attached to glycerol backbone"

        if ester_count < 2:
            return False, f"Found {ester_count} ester groups connected to glycerol carbons, need at least 2"

        # Check passed
        return True, "Contains glycerol backbone with phosphate group and two fatty acid chains attached via ester bonds"

    return False, "Does not match phosphatidic acid structure"

__metadata__ = {   'chemical_class': {   'id': 'CHEBI:49105',
                              'name': 'phosphatidic acid',
                              'definition': 'A derivative of glycerol in which one hydroxy group, commonly but not necessarily primary, is esterified with phosphoric acid and the other two are esterified with fatty acids.',
                              'parents': []},
        'config': {   'llm_model_name': 'lbl/claude-sonnet',
                      'f1_threshold': 0.8,
                      'max_attempts': 5},
        'message': None,
        'attempt': 1,
        'success': True,
        'best': True,
        'error': '',
        'stdout': None}