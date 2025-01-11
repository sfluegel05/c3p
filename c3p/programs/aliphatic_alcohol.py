"""
Classifies: CHEBI:2571 aliphatic alcohol
"""
"""
Classifies: aliphatic alcohol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_aliphatic_alcohol(smiles: str):
    """
    Determines if a molecule is an aliphatic alcohol based on its SMILES string.
    An aliphatic alcohol is an alcohol derived from an aliphatic compound.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aliphatic alcohol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for hydroxyl groups (excluding carboxylic acids)
    alcohol_pattern = Chem.MolFromSmarts("[CX4][OH1]")
    if not mol.HasSubstructMatch(alcohol_pattern):
        return False, "No aliphatic hydroxyl groups found"
    
    # Count hydroxyl groups
    hydroxyl_matches = len(mol.GetSubstructMatches(alcohol_pattern))
    
    # Count carboxylic acid groups
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=O)[OH1]")
    carboxyl_matches = len(mol.GetSubstructMatches(carboxyl_pattern))
    
    # If all OH groups are from carboxylic acids, it's not an alcohol
    if hydroxyl_matches == 0:
        return False, "No alcohol groups present"
    
    # Calculate aromaticity
    aromatic_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
    total_atoms = mol.GetNumAtoms()
    aromatic_ratio = aromatic_atoms / total_atoms if total_atoms > 0 else 0
    
    # If the molecule is primarily aromatic (>50%), reject it
    if aromatic_ratio > 0.5:
        return False, "Molecule is primarily aromatic"
    
    # Get all carbons attached to OH groups
    carbons_with_oh = []
    for match in mol.GetSubstructMatches(alcohol_pattern):
        c_idx = match[0]  # Carbon index
        carbons_with_oh.append(mol.GetAtomWithIdx(c_idx))
    
    # Check if any carbon with OH is sp3 hybridized
    has_sp3_alcohol = False
    for carbon in carbons_with_oh:
        if carbon.GetHybridization() == Chem.HybridizationType.SP3:
            has_sp3_alcohol = True
            break
            
    if not has_sp3_alcohol:
        return False, "No sp3 carbons with hydroxyl groups found"
    
    # Success case
    if hydroxyl_matches == 1:
        return True, "Contains one aliphatic hydroxyl group"
    else:
        return True, f"Contains {hydroxyl_matches} aliphatic hydroxyl groups"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:2571',
                          'name': 'aliphatic alcohol',
                          'definition': 'An  alcohol derived from an aliphatic '
                                        'compound.',
                          'parents': ['CHEBI:30879'],
                          'xrefs': ['KEGG:C02525'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'c1c[nH+]c[nH]1',
                                     'name': 'imidazolium cation',
                                     'reason': 'No aliphatic hydroxyl groups '
                                               'found'},
                                 {   'smiles': 'ClC(Cl)[C@H](O)CC=1OC(=O)C=2C(O)=CC(=CC2C1)O',
                                     'name': 'Desmethyldichlorodiaportin',
                                     'reason': 'Molecule is primarily '
                                               'aromatic'},
                                 {   'smiles': 'ClC=1C(=C(O)C2=C(C1O)C(=O)C=CC2=O)CC=C(C)C',
                                     'name': 'Chlorosesamone',
                                     'reason': 'No aliphatic hydroxyl groups '
                                               'found'},
                                 {   'smiles': 'CCCCC[C@H]1O[C@H]1C\\C=C/CCCCCCCC(O)=O',
                                     'name': '(+)-vernolic acid',
                                     'reason': 'No aliphatic hydroxyl groups '
                                               'found'},
                                 {   'smiles': 'O(C1=C(C=CC=C1C(O)=O)C)C(=O)C',
                                     'name': 'CRESOPYRINE',
                                     'reason': 'No aliphatic hydroxyl groups '
                                               'found'},
                                 {   'smiles': 'O1C2C3C(CCC3=C)C(CCC2C(C1=O)=C)=C',
                                     'name': '3,6,9-Trimethylidene-3a,4,5,6a,7,8,9a,9b-octahydroazuleno[4,5-b]furan-2-one',
                                     'reason': 'No aliphatic hydroxyl groups '
                                               'found'},
                                 {   'smiles': 'S(C(S)CC)CCC',
                                     'name': 'xi-1-(Propylthio)-1-propanethiol',
                                     'reason': 'No aliphatic hydroxyl groups '
                                               'found'},
                                 {   'smiles': 'O=C(OC)[C@]1(C=C(CCC1)CC)CC(=O)O',
                                     'name': '(-)-petchioic B',
                                     'reason': 'No aliphatic hydroxyl groups '
                                               'found'},
                                 {   'smiles': 'C[C@@H]1CN(C(=O)C2=C(N=CC(=C2)C3=CC=CC=C3C)O[C@@H]1CN(C)CC4=CN=CN=C4)[C@@H](C)CO',
                                     'name': '(2S,3R)-5-[(2S)-1-hydroxypropan-2-yl]-3-methyl-8-(2-methylphenyl)-2-[[methyl(5-pyrimidinylmethyl)amino]methyl]-3,4-dihydro-2H-pyrido[2,3-b][1,5]oxazocin-6-one',
                                     'reason': 'Molecule is primarily '
                                               'aromatic'},
                                 {   'smiles': 'O=C(N[C@@H](CC=1C=2C(NC1)=CC=CC2)C(O)=O)[C@@H](NC(=O)[C@@H](N)CCCN=C(N)N)CC(C)C',
                                     'name': 'Arg-Leu-Trp',
                                     'reason': 'No aliphatic hydroxyl groups '
                                               'found'}],
    'sample_false_negatives': [],
    'attempt': 0,
    'reasoning': 'Let me create a program to classify aliphatic alcohols. '
                 "Here's my reasoning:\n"
                 '\n'
                 '1. An aliphatic alcohol should have:\n'
                 '   - At least one hydroxyl (-OH) group\n'
                 '   - A primarily aliphatic (non-aromatic) structure\n'
                 '   - The OH group should be attached to an sp3 carbon (not '
                 'part of carboxylic acid)\n'
                 '\n'
                 '2. We should exclude:\n'
                 '   - Molecules without OH groups\n'
                 '   - Molecules where OH is only part of carboxylic acids\n'
                 '   - Molecules where OH is primarily attached to aromatic '
                 'systems\n'
                 '\n'
                 '3. Key implementation points:\n'
                 '   - Check for presence of OH groups\n'
                 '   - Check if the carbon attached to OH is sp3 hybridized\n'
                 '   - Calculate aromaticity to ensure the molecule is '
                 'primarily aliphatic\n'
                 '   - Handle edge cases and invalid SMILES\n'
                 '\n'
                 "Here's the implementation:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 731,
    'num_false_positives': 100,
    'num_true_negatives': 118,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.8796630565583634,
    'recall': 1.0,
    'f1': 0.9359795134443022,
    'accuracy': 0.8946259220231823,
    'negative_predictive_value': 1.0}