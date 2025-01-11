"""
Classifies: CHEBI:64611 ether lipid
"""
"""
Classifies: CHEBI:12307 ether lipid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_ether_lipid(smiles: str):
    """
    Determines if a molecule is an ether lipid based on its SMILES string.
    An ether lipid is similar in structure to a glycerolipid but in which one or more of the carbon atoms on glycerol
    is bonded to an alkyl chain via an ether linkage, as opposed to the usual ester linkage.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an ether lipid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone pattern (three carbons connected in sequence)
    glycerol_pattern = Chem.MolFromSmarts("C-C-C")
    matches = mol.GetSubstructMatches(glycerol_pattern)
    if not matches:
        return False, "No glycerol backbone found"

    # For each glycerol backbone found
    for match in matches:
        c1_idx, c2_idx, c3_idx = match
        glycerol_carbon_indices = set([c1_idx, c2_idx, c3_idx])

        # Check for ether linkage on glycerol carbons
        ether_found = False
        for c_idx in glycerol_carbon_indices:
            atom = mol.GetAtomWithIdx(c_idx)
            # Iterate over neighbors to find oxygen connected via single bond
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8:
                    bond = mol.GetBondBetweenAtoms(c_idx, neighbor.GetIdx())
                    if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                        oxygen_atom = neighbor
                        # Check if oxygen is connected to another carbon not in glycerol backbone
                        for o_neighbor in oxygen_atom.GetNeighbors():
                            if o_neighbor.GetIdx() != c_idx and o_neighbor.GetAtomicNum() == 6:
                                if o_neighbor.GetIdx() not in glycerol_carbon_indices:
                                    o_bond = mol.GetBondBetweenAtoms(oxygen_atom.GetIdx(), o_neighbor.GetIdx())
                                    if o_bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                                        # Found ether linkage
                                        ether_found = True
                                        break
                        if ether_found:
                            break
            if ether_found:
                break

        if ether_found:
            return True, "Contains glycerol backbone with ether linkage"

    return False, "No ether linkage found on glycerol backbone"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:64611',
                          'name': 'ether lipid',
                          'definition': 'A lipid similar in structure to a '
                                        'glycerolipid but in which one or more '
                                        'of the carbon atoms on glycerol is '
                                        'bonded to an alkyl chain via an ether '
                                        'linkage, as opposed to the usual '
                                        'ester linkage.',
                          'parents': ['CHEBI:18059', 'CHEBI:52575'],
                          'xrefs': [   'PMID:21309516',
                                       'PMID:22114698',
                                       'PMID:22148427',
                                       'PMID:22306069',
                                       'PMID:22348073',
                                       'PMID:22366205',
                                       'PMID:22506086',
                                       'PMID:22609598',
                                       'Wikipedia:Ether_lipid'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'CNC(=O)CC[C@H](N)C(O)=O',
                                     'name': 'N(5)-methyl-L-glutamine',
                                     'reason': 'No ether linkage found on '
                                               'glycerol backbone'},
                                 {   'smiles': '[H][C@@]1(CC[C@]2(C)[C@]1([H])CC[C@]1([H])[C@@]3(C)CCCC(C)(C)[C@]3([H])CC[C@@]21C)C(=C)CCC=C(C)C',
                                     'name': 'dammara-20,24-diene',
                                     'reason': 'No ether linkage found on '
                                               'glycerol backbone'},
                                 {   'smiles': 'O=C1C2=C(O)C3=C(O)C=CC=C3C(=C2C(C(=O)C)C(C1)(O)C)C=4C=5C(C(O)=C6C4C(C(=O)C)C(O)(C)CC6=O)=C(O)C=CC5',
                                     'name': 'A 39183A',
                                     'reason': 'No ether linkage found on '
                                               'glycerol backbone'},
                                 {   'smiles': 'C[C@H](O)[C@@H](O)C1=Nc2c(NC1)[nH]c(N)nc2=O',
                                     'name': 'L-threo-7,8-dihydrobiopterin',
                                     'reason': 'No ether linkage found on '
                                               'glycerol backbone'},
                                 {   'smiles': 'C(CCN1CCCCC1)(CCN(C(C)C)C(C)=O)(C(N)=O)C2=C(Cl)C=CC=C2',
                                     'name': 'bidisomide',
                                     'reason': 'No ether linkage found on '
                                               'glycerol backbone'},
                                 {   'smiles': 'OCCCCCCCCCCC#CC=C',
                                     'name': '13-Tetradece-11-yn-1-ol',
                                     'reason': 'No ether linkage found on '
                                               'glycerol backbone'},
                                 {   'smiles': 'CC1=CC=CC=C1C2=CC=C(C=C2)C3[C@H]4CNC[C@@H]3N4C(=O)NC5=CC=C(C=C5)Cl',
                                     'name': '(1S,5R)-N-(4-chlorophenyl)-7-[4-(2-methylphenyl)phenyl]-3,6-diazabicyclo[3.1.1]heptane-6-carboxamide',
                                     'reason': 'No ether linkage found on '
                                               'glycerol backbone'},
                                 {   'smiles': 'OC(=O)c1cccc(Oc2ccccc2)c1',
                                     'name': '3-phenoxybenzoic acid',
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': 'N1(C(C(=O)C=2C1=CC=C(C#N)C2)=O)CC3=C(ON=C3C)C',
                                     'name': '1-[(3,5-Dimethyl-1,2-oxazol-4-yl)methyl]-2,3-dioxoindole-5-carbonitrile',
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': 'CCCN',
                                     'name': 'propylamine',
                                     'reason': 'No ether linkage found on '
                                               'glycerol backbone'}],
    'sample_false_negatives': [   {   'smiles': 'O[C@@H]1C/C(=C\\C=C\\2/[C@]3([C@@]([C@](CC3)([C@H](C)C=C[C@@H](C(C)C)C)[H])(CCC2)C)[H])/C([C@@H](O)C1)=C',
                                      'name': '1-octadecyl lysophosphatidic '
                                              'acid',
                                      'reason': 'No ether linkage found on '
                                                'glycerol backbone'}],
    'attempt': 0,
    'reasoning': 'To determine if a molecule is an ether lipid, the code '
                 'will:\n'
                 '\n'
                 '- Parse the SMILES string into a molecule.\n'
                 '- Identify the glycerol backbone, which is a sequence of '
                 'three connected carbon atoms (C-C-C).\n'
                 '- For each carbon atom in the glycerol backbone, check if it '
                 'is connected via a single bond to an oxygen atom.\n'
                 '- Check if that oxygen atom is connected via a single bond '
                 'to a carbon atom not part of the glycerol backbone '
                 '(indicating an ether linkage to an alkyl chain).\n'
                 '- If at least one glycerol carbon has such an ether linkage, '
                 'classify the molecule as an ether lipid.\n'
                 '\n'
                 'Here is the code:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 251,
    'num_false_positives': 100,
    'num_true_negatives': 113,
    'num_false_negatives': 1,
    'num_negatives': None,
    'precision': 0.7150997150997151,
    'recall': 0.996031746031746,
    'f1': 0.8325041459369819,
    'accuracy': 0.7827956989247312,
    'negative_predictive_value': 0.9912280701754386}