"""
Classifies: CHEBI:26267 proanthocyanidin
"""
"""
Classifies: CHEBI:38142 proanthocyanidin

A proanthocyanidin is a flavonoid oligomer obtained by the the condensation of two or more units of hydroxyflavans.
Common features:
- Contains flavan-3-ol units (e.g., catechin, epicatechin, gallocatechin, epigallocatechin)
- Units linked via C-C or ether bonds, often 4→8, 4→6, or 2→O→7 linkages
- May contain galloyl or other ester groups
- Typically higher molecular weight oligomers or polymers
"""

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_proanthocyanidin(smiles: str):
    """
    Determines if a molecule is a proanthocyanidin based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a proanthocyanidin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for flavan-3-ol units
    flavan_pattern = Chem.MolFromSmarts("[C&r5]1=C([C@@H](O)[C@H](O)[C@@]2(O)Cc3c(O)cc(O)cc3O2)c2ccccc12")
    flavan_matches = mol.GetSubstructMatches(flavan_pattern)
    if not flavan_matches:
        return False, "No flavan-3-ol units found"
    
    # Look for common linkages (4→8, 4→6, 2→O→7)
    link_patterns = [
        Chem.MolFromSmarts("[C&r4]@[C&r8]"),  # 4→8 linkage
        Chem.MolFromSmarts("[C&r4]@[C&r6]"),  # 4→6 linkage
        Chem.MolFromSmarts("[C&r2]@[O&r7]")   # 2→O→7 linkage
    ]
    link_matches = []
    for pattern in link_patterns:
        link_matches.extend(mol.GetSubstructMatches(pattern))
    if not link_matches:
        return False, "No characteristic linkages found between flavan units"
    
    # Check for galloyl or other ester groups
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    # Check molecular weight - proanthocyanidins typically >500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        if not ester_matches:
            return False, "Molecular weight too low for proanthocyanidin"
    
    return True, "Contains flavan-3-ol units linked via characteristic bonds, potentially with ester groups"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:38142',
        'name': 'proanthocyanidin',
        'definition': 'A flavonoid oligomer obtained by the the condensation of two or more units of hydroxyflavans.',
        'parents': ['CHEBI:24034', 'CHEBI:27856']
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
    'num_true_positives': 130,
    'num_false_positives': 8,
    'num_true_negatives': 182359,
    'num_false_negatives': 1,
    'num_negatives': None,
    'precision': 0.9418604651162791,
    'recall': 0.9923076923076923,
    'f1': 0.9667769727484958,
    'accuracy': 0.9995654976664447
}