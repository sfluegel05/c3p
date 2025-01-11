"""
Classifies: CHEBI:32863 secondary amine
"""
"""
Classifies: CHEBI:secondary_amine
"""
from rdkit import Chem

def is_secondary_amine(smiles: str):
    """
    Determines if a molecule is a secondary amine based on its SMILES string.
    A secondary amine is characterized by having a nitrogen atom bonded to two hydrocarbyl groups and one hydrogen.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary amine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Iterate over all atoms and find secondary amine nitrogens
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:  # Check for nitrogen atoms
            # Count the number of carbon bonds
            carbon_bond_count = sum(1 for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() == 6)
            non_hydrocarbon_bond_count = sum(1 for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() != 6)
            
            # Check if there are exactly two carbon bonds
            if carbon_bond_count == 2 and non_hydrocarbon_bond_count <= 1:
                return True, "Contains nitrogen with two hydrocarbyl groups"

    return False, "No nitrogen atom with two hydrocarbyl groups found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:32863',
                          'name': 'secondary amine',
                          'definition': 'A compound formally derived from '
                                        'ammonia by replacing two hydrogen '
                                        'atoms by hydrocarbyl groups.',
                          'parents': ['CHEBI:32952', 'CHEBI:50995'],
                          'xrefs': ['KEGG:C02324'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'O=C1C2=C(OC(=C1)C)C3=C(OC)C=C(OC)C=C3C(=C2O)C4=C5OC(=CC(C5=C(O)C=6C4=CC(OC)=CC6OC)=O)C',
                                     'name': 'Isonigerone',
                                     'reason': 'No nitrogen atom with two '
                                               'hydrocarbyl groups found'},
                                 {   'smiles': 'CCCCCCCCCCCCCC(=O)OC[C@H](COP([O-])(=O)OC[C@@H](O)CO)OC(=O)CCCCCCC\\C=C/CCCCCCCC',
                                     'name': '1-myristoyl-2-oleoyl-sn-glycero-3-phosphatidylglycerol(1-)',
                                     'reason': 'No nitrogen atom with two '
                                               'hydrocarbyl groups found'},
                                 {   'smiles': 'S([C@H]1N(C(=O)/C(=C\\C2=CC=C(OCC=C(C)C)C=C2)/N(C1=O)C)C)C',
                                     'name': 'Fusaperazine F',
                                     'reason': 'No nitrogen atom with two '
                                               'hydrocarbyl groups found'},
                                 {   'smiles': 'Oc1ccccc1I',
                                     'name': '2-iodophenol',
                                     'reason': 'No nitrogen atom with two '
                                               'hydrocarbyl groups found'},
                                 {   'smiles': 'O(C1=C(OC)C=C(C2=C(O)C(OC)=C(C3=CC=CC=C3)C=C2OC)C=C1)C/C=C(/CO)\\C',
                                     'name': 'Prenylterphenyllin F',
                                     'reason': 'No nitrogen atom with two '
                                               'hydrocarbyl groups found'},
                                 {   'smiles': 'O1[C@@H]([C@@H](O)[C@H](O[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2O)CO)[C@@H](O)[C@@H]1OC[C@H]3O[C@@H](OC[C@H]4O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]4O)[C@H](O)[C@@H](O)[C@@H]3O)CO[C@@H]5O[C@@H]([C@@H](O)[C@H](O[C@@H]6O[C@@H]([C@@H](O)[C@H](O)[C@H]6O)CO)[C@H]5O)CO[C@@H]7O[C@@H]([C@@H](O)[C@H](O)[C@H]7O)CO',
                                     'name': '(2R,3R,4S,5S,6R)-6-[[(2R,3R,4S,5S,6R)-6-[[(2R,3R,4S,5R,6R)-6-[[(2R,3R,4S,5R,6R)-3,5-Dihydroxy-4-[(2R,3R,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-6-[[(2R,3R,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxymethyl]oxan-2-yl]oxymethyl]-3,5-dihydroxy-4-[(2R,3R,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxyoxan-2-yl]oxymethyl]-3,4,5-trihydroxyoxan-2-yl]oxymethyl]oxane-2,3,4,5-tetrol',
                                     'reason': 'No nitrogen atom with two '
                                               'hydrocarbyl groups found'},
                                 {   'smiles': 'P(OC[C@H](OC(=O)CCCCCCCCCCCCCCC)COC(=O)CCCCCCCCCCCCCCCC)(OC[C@H](N)C(O)=O)(O)=O',
                                     'name': 'PS(17:0/16:0)',
                                     'reason': 'No nitrogen atom with two '
                                               'hydrocarbyl groups found'},
                                 {   'smiles': 'C[C@@H]([C@H]1CC[C@H]2[C@@H]3CC=C4C=CC(=O)[C@]4(C)[C@H]3CC[C@]12C)[C@H]1CC(C)=C(CO)C(=O)O1',
                                     'name': 'withalongolide F',
                                     'reason': 'No nitrogen atom with two '
                                               'hydrocarbyl groups found'},
                                 {   'smiles': 'OC[C@@H](O)[C@@H](O)[C@H](O)[C@H](O)CC(=O)C([O-])=O',
                                     'name': 'keto-3-deoxy-D-manno-octulosonate',
                                     'reason': 'No nitrogen atom with two '
                                               'hydrocarbyl groups found'},
                                 {   'smiles': 'O=C1C(=C(OC(=C1C)[C@@H](O)CCCCCC)OC)C',
                                     'name': 'Nocapyrone J',
                                     'reason': 'No nitrogen atom with two '
                                               'hydrocarbyl groups found'}],
    'sample_false_negatives': [],
    'attempt': 0,
    'reasoning': 'To classify chemical compounds as secondary amines based on '
                 'their SMILES strings, we should focus on the definition of '
                 'secondary amines: compounds derived from ammonia by '
                 'replacing two hydrogen atoms with hydrocarbyl groups. This '
                 'means that we are looking for compounds containing a '
                 'nitrogen atom bonded to exactly two carbon atoms.\n'
                 '\n'
                 "Here's a step-by-step breakdown to achieve this:\n"
                 '\n'
                 '1. **Parse the SMILES String**: Utilize RDKit to convert the '
                 'SMILES string into a molecular object.\n'
                 '2. **Identify Nitrogen Atoms**: Locate all nitrogen atoms in '
                 'the molecule.\n'
                 '3. **Determine Connectivity**: For each nitrogen atom, check '
                 'if it is bonded to two carbon atoms and possibly a third '
                 'atom (like hydrogen or something else).\n'
                 '4. **Classification**: If 1+ nitrogen atoms in the molecule '
                 'meet the criteria for secondary amines, classify '
                 'accordingly.\n'
                 '\n'
                 'Below is the implementation using RDKit.',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 44,
    'num_false_positives': 100,
    'num_true_negatives': 87,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.3055555555555556,
    'recall': 1.0,
    'f1': 0.46808510638297873,
    'accuracy': 0.5670995670995671,
    'negative_predictive_value': 1.0}