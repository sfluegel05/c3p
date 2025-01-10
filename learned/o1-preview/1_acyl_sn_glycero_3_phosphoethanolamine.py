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
    esterified at the sn-1 position, a hydroxyl group at the sn-2 position, and a phosphoethanolamine group at the sn-3 position,
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
    # Glycerol backbone with ester at sn-1, hydroxyl at sn-2, phosphoethanolamine at sn-3
    glycerol_pattern = Chem.MolFromSmarts("""
        [C@@H]([O])[C@@H](COP(=O)(OCCN)O)O
        """)
    ester_pattern = Chem.MolFromSmarts("""
        [C](=O)O[C@@H]
        """)
    phospho_ethanolamine_pattern = Chem.MolFromSmarts("""
        P(=O)(O[C][C]N)O
        """)

    # Search for ester linkage at sn-1 position
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester linkage found at sn-1 position"

    # Search for glycerol backbone with hydroxyl at sn-2 and phosphoethanolamine at sn-3
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)
    if not glycerol_matches:
        return False, "No glycerol backbone with required functional groups found"

    # Search for phosphoethanolamine group
    phospho_matches = mol.GetSubstructMatches(phospho_ethanolamine_pattern)
    if not phospho_matches:
        return False, "No phosphoethanolamine group found at sn-3 position"

    # Optional: Check chiral configuration at sn-2 (C2)
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    sn2_chiral = None
    for idx, config in chiral_centers:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetSymbol() == 'C':
            # Check if this carbon is connected to OH and part of glycerol backbone
            neighbors = [n.GetSymbol() for n in atom.GetNeighbors()]
            if neighbors.count('O') >= 1 and neighbors.count('C') >= 2:
                sn2_chiral = config
                break
    if sn2_chiral:
        if sn2_chiral != 'R':
            return False, f"Chiral center at sn-2 position is not (R)-configuration (found {sn2_chiral})"

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
    'attempt': 2,
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