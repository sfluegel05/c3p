"""
Classifies: CHEBI:68452 azole
"""
"""
Classifies: azole compounds
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_azole(smiles: str):
    """
    Determines if a molecule contains an azole ring based on its SMILES string.
    An azole is a five-membered aromatic heterocycle containing at least one nitrogen atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains an azole ring, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate 2D coordinates for ring detection
    AllChem.Compute2DCoords(mol)

    # Find all rings in the molecule
    rings = mol.GetRingInfo()
    
    # Check each ring of size 5
    for ring_atoms in rings.AtomRings():
        if len(ring_atoms) != 5:
            continue
            
        # Get the atoms in the ring
        ring_atom_objects = [mol.GetAtomWithIdx(i) for i in ring_atoms]
        
        # Check if ring contains at least one nitrogen
        has_nitrogen = False
        other_heteroatoms = False
        all_aromatic = True
        
        for atom in ring_atom_objects:
            # Check for nitrogen
            if atom.GetAtomicNum() == 7:
                has_nitrogen = True
            # Check for other heteroatoms
            elif atom.GetAtomicNum() in [8, 16]:  # O or S
                other_heteroatoms = True
            # Check if all atoms are aromatic
            if not atom.GetIsAromatic():
                all_aromatic = False
                
        # If we found a five-membered aromatic ring with at least one nitrogen
        if has_nitrogen and all_aromatic:
            msg = "Contains a five-membered aromatic ring with nitrogen"
            if other_heteroatoms:
                msg += " and additional heteroatoms"
            return True, msg

    return False, "No azole ring found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:68452',
                          'name': 'azole',
                          'definition': 'Any monocyclic heteroarene consisting '
                                        'of a five-membered ring containing '
                                        'nitrogen. Azoles can also contain one '
                                        'or more other non-carbon atoms, such '
                                        'as nitrogen, sulfur or oxygen.',
                          'parents': ['CHEBI:38101', 'CHEBI:38179'],
                          'xrefs': ['Wikipedia:Azole'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [],
    'sample_false_negatives': [   {   'smiles': 'C1CN=C(Cc2ccccc2)N1',
                                      'name': 'tolazoline',
                                      'reason': 'No azole ring found'},
                                  {   'smiles': 'C[C@@](O)([C@H]1[C@H](O)C[C@@]2(C)[C@@H]3CC=C4[C@@H](C=C(O[C@@H]5O[C@H](CO)[C@@H](O)[C@H](O)[C@H]5O)C(=O)C4(C)C)[C@]3(C)C(=O)C[C@]12C)C(=O)C[C@H]1NC(=O)OC1(C)C',
                                      'name': '(16alpha,20R,24R)-24N,25-carbamoyloxy-2,16,20-trihydroxycucurbita-1,5-diene-3,11,22-trione '
                                              '2-O-beta-D-glucopyranoside',
                                      'reason': 'No azole ring found'},
                                  {   'smiles': 'COC1=C(Br)[C@H](O)[C@]2(CC(=NO2)C(=O)NCCCCNC(=O)C2=NO[C@@]3(C2)C=C(Br)C(OC)=C(Br)[C@@H]3O)C=C1Br',
                                      'name': 'aerothionin',
                                      'reason': 'No azole ring found'},
                                  {   'smiles': 'CCCCCCCC(=O)C(C)C(=O)N1CCC=C1',
                                      'name': '1-(2,3-dihydro-1H-pyrrol-1-yl)-2-methyldecane-1,3-dione',
                                      'reason': 'No azole ring found'},
                                  {   'smiles': '[H][C@]1(N=C(O[C@H]1C)c1cccc(O)c1O)C(=O)NCCCN(CCCNC(=O)c1cccc(O)c1O)C(=O)[C@@]1([H])N=C(O[C@@H]1C)c1cccc(O)c1O',
                                      'name': 'vibriobactin',
                                      'reason': 'No azole ring found'},
                                  {   'smiles': 'CO[C@H](CC\\C(C)=C\\C=C\\CC\\C=C/[C@@H]1CSC(=N1)[C@@H]1C[C@@H]1C)CC=C',
                                      'name': 'Curacin A',
                                      'reason': 'No azole ring found'},
                                  {   'smiles': 'CC(C1=NCCS1)=O',
                                      'name': '2-Acetyl-4,5-dihydrothiazole',
                                      'reason': 'No azole ring found'},
                                  {   'smiles': '[H]C(OCC)=C1N=C(OC1=O)c1ccco1',
                                      'name': '4-(ethoxymethylene)-2-(2-furyl)oxazol-5-one',
                                      'reason': 'No azole ring found'},
                                  {   'smiles': 'CCCCCCCCCCCCCCCCC\\C=C/C(=O)NCCCCC(NC(=O)C1(C)COC(=N1)c1ccccc1O)C(=O)OC(C)CC(=O)NC1CCCCNC1=O',
                                      'name': 'DDM-838',
                                      'reason': 'No azole ring found'},
                                  {   'smiles': 'C=1(C(=CC(=CN1)C)C(O)=O)C=2NC([C@@](N2)(C)C(C)C)=O',
                                      'name': '(R)-imazapic',
                                      'reason': 'No azole ring found'}],
    'attempt': 0,
    'reasoning': 'Let me create a program to classify azoles. The key '
                 'characteristics are:\n'
                 '1. Five-membered heterocyclic ring\n'
                 '2. Must contain at least one nitrogen\n'
                 '3. May contain additional heteroatoms (N, S, O)\n'
                 '4. Must be aromatic\n'
                 '\n'
                 "Here's the implementation:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 1628,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 100,
    'num_negatives': None,
    'precision': 1.0,
    'recall': 0.9421296296296297,
    'f1': 0.970202622169249,
    'accuracy': 0.9421296296296297,
    'negative_predictive_value': 0.0}