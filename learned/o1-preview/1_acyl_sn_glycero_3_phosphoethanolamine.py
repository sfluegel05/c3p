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
    This class has a glycerol backbone with an acyl group esterified at the sn-1 position,
    a hydroxyl group at the sn-2 position, and a phosphoethanolamine group at the sn-3 position,
    with (R)-configuration at the sn-2 carbon if stereochemistry is specified.

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

    # Assign stereochemistry
    Chem.AssignAtomChiralTagsFromStructure(mol)

    # Define SMARTS patterns
    # Ester linkage at sn-1 position
    ester_sn1_pattern = Chem.MolFromSmarts("[$([O][C](=O)[C])]")

    # Hydroxyl group at sn-2 position
    hydroxyl_sn2_pattern = Chem.MolFromSmarts("[C@H](O)")

    # Phosphoethanolamine group at sn-3 position
    phospho_ethanolamine_pattern = Chem.MolFromSmarts("[O][P](=O)([O])[O][C][C][N]")

    # Validate patterns
    if None in (ester_sn1_pattern, hydroxyl_sn2_pattern, phospho_ethanolamine_pattern):
        return False, "Error in SMARTS pattern definitions"

    # Search for ester linkage at sn-1 position
    ester_matches = mol.GetSubstructMatches(ester_sn1_pattern)
    if not ester_matches:
        return False, "No ester linkage found at sn-1 position"

    # Search for hydroxyl group at sn-2 position
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_sn2_pattern)
    if not hydroxyl_matches:
        return False, "No hydroxyl group found at sn-2 position"

    # Search for phosphoethanolamine group at sn-3 position
    phospho_matches = mol.GetSubstructMatches(phospho_ethanolamine_pattern)
    if not phospho_matches:
        return False, "No phosphoethanolamine group found at sn-3 position"

    # Identify sn-2 chiral carbon
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=False)
    sn2_chiral = None
    for idx, config in chiral_centers:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetSymbol() == 'C' and config == 'R':
            # Check if this carbon matches the hydroxyl pattern
            if mol.HasSubstructMatch(Chem.MolFromSmarts("[C@H](O)")):
                sn2_chiral = idx
                break

    if sn2_chiral is None:
        return False, "No chiral center with (R)-configuration at sn-2 position"

    # Verify connectivity of sn-2 carbon
    sn2_atom = mol.GetAtomWithIdx(sn2_chiral)
    neighbor_elements = [nbr.GetSymbol() for nbr in sn2_atom.GetNeighbors()]
    if neighbor_elements.count('O') < 1 or neighbor_elements.count('C') < 3:
        return False, "sn-2 carbon does not have the correct connectivity"

    return True, "Molecule is a 1-acyl-sn-glycero-3-phosphoethanolamine"

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
    'attempt': 3,
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