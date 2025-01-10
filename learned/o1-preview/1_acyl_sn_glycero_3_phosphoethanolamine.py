"""
Classifies: CHEBI:29017 1-acyl-sn-glycero-3-phosphoethanolamine
"""
"""
Classifies: 1-acyl-sn-glycero-3-phosphoethanolamine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_1_acyl_sn_glycero_3_phosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-acyl-sn-glycero-3-phosphoethanolamine based on its SMILES string.
    A 1-acyl-sn-glycero-3-phosphoethanolamine has a glycerol backbone with an acyl group
    esterified at the sn-1 position, a hydroxyl group at sn-2 position, and a phosphoethanolamine group at the sn-3 position,
    with (R)-configuration at the sn-2 position if stereochemistry is specified.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-acyl-sn-glycero-3-phosphoethanolamine, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add hydrogens to correctly perceive chiral centers
    mol = Chem.AddHs(mol)

    # Define SMARTS patterns
    glycerol_pattern = Chem.MolFromSmarts("[C;!R]-[C;!R]-[C;!R]")
    ester_pattern = Chem.MolFromSmarts("C(=O)[O][C]")
    hydroxyl_pattern = Chem.MolFromSmarts("[C][O;H1]")
    phospho_ethanolamine_pattern = Chem.MolFromSmarts("P(=O)([O-])[O][C][C][N]")

    # Find glycerol backbone
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)
    if not glycerol_matches:
        return False, "No glycerol backbone found"

    for match in glycerol_matches:
        c1_idx, c2_idx, c3_idx = match
        c1 = mol.GetAtomWithIdx(c1_idx)
        c2 = mol.GetAtomWithIdx(c2_idx)
        c3 = mol.GetAtomWithIdx(c3_idx)

        # Check for ester linkage at sn-1 (C1)
        ester_found = False
        for neighbor in c1.GetNeighbors():
            if neighbor.GetIdx() != c2_idx:
                bond = mol.GetBondBetweenAtoms(c1.GetIdx(), neighbor.GetIdx())
                if bond.GetBondType() == Chem.rdchem.BondType.SINGLE and neighbor.GetAtomicNum() == 8:
                    # Check if this oxygen is part of an ester linkage
                    for ester_match in mol.GetSubstructMatches(ester_pattern):
                        if neighbor.GetIdx() == ester_match[2] and c1.GetIdx() == ester_match[3]:
                            ester_found = True
                            break
            if ester_found:
                break
        if not ester_found:
            continue

        # Check for hydroxyl group at sn-2 (C2)
        hydroxyl_found = False
        for neighbor in c2.GetNeighbors():
            if neighbor.GetIdx() != c1_idx and neighbor.GetIdx() != c3_idx:
                if mol.GetBondBetweenAtoms(c2.GetIdx(), neighbor.GetIdx()).GetBondType() == Chem.rdchem.BondType.SINGLE and neighbor.GetAtomicNum() == 8:
                    if neighbor.GetTotalNumHs() > 0:
                        hydroxyl_found = True
                        break
        if not hydroxyl_found:
            continue

        # Check for phosphoethanolamine group at sn-3 (C3)
        phospho_found = False
        for neighbor in c3.GetNeighbors():
            if neighbor.GetIdx() != c2_idx and neighbor.GetAtomicNum() == 8:
                # Check if this oxygen is connected to phosphoethanolamine
                for phospho_match in mol.GetSubstructMatches(phospho_ethanolamine_pattern):
                    if neighbor.GetIdx() == phospho_match[2]:
                        phospho_found = True
                        break
            if phospho_found:
                break
        if not phospho_found:
            continue

        # Optional: Check chiral configuration at sn-2 (C2)
        chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
        sn2_chiral = None
        for idx, config in chiral_centers:
            if idx == c2_idx:
                sn2_chiral = config
                break
        if sn2_chiral:
            if sn2_chiral != 'R':
                return False, f"Chiral center at sn-2 position is not (R)-configuration (found {sn2_chiral})"

        # All checks passed
        return True, "Molecule is a 1-acyl-sn-glycero-3-phosphoethanolamine"

    return False, "Glycerol backbone with required functional groups not found"

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': '1-acyl-sn-glycero-3-phosphoethanolamine',
        'definition': 'A 1-O-acylglycerophosphoethanolamine having (R)-configuration.',
        'parents': []
    },
    'config': {
        'llm_model_name': 'lbl/claude-sonnet',
        'f1_threshold': 0.8,
        'max_attempts': 5,
        'max_positive_instances': None,
        'max_positive_to_test': None,
        'max_negative_to_test': None,
        'max_positive_in_prompt': 50,
        'max_negative_in_prompt': 20,
        'max_instances_in_prompt': 100,
        'test_proportion': 0.1
    },
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
    'accuracy': None
}