"""
Classifies: CHEBI:26004 phenylpropanoid
"""
"""
Classifies: phenylpropanoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_phenylpropanoid(smiles: str):
    """
    Determines if a molecule is a phenylpropanoid based on its SMILES string.
    A phenylpropanoid is any organic aromatic compound with a structure based on a phenylpropane skeleton.
    This includes compounds like flavonoids, coumarins, lignins, stilbenes, and others with a C6-C3 backbone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phenylpropanoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define generalized SMARTS patterns for phenylpropanoid structures
    patterns = {
        'phenylpropane_core': Chem.MolFromSmarts('c1ccccc1CCC'),  # Phenyl ring linked to 3-carbon chain
        'phenylpropene_core': Chem.MolFromSmarts('c1ccccc1C=CC'),  # Phenyl ring linked to propenyl chain
        'phenylpropyne_core': Chem.MolFromSmarts('c1ccccc1C#CC'),  # Phenyl ring linked to propynyl chain
        'cinnamic_acid': Chem.MolFromSmarts('c1ccccc1C=CC(=O)O'),  # Cinnamic acid scaffold
        'coumarin_core': Chem.MolFromSmarts('O=C1C=CC2=CC=CC=C2O1'),  # Coumarin core
        'flavonoid_core': Chem.MolFromSmarts('c1cc(c(cc1)-c1coc2c1ccc(=O)c(=O)c2)O'),  # Flavonoid skeleton
        'isoflavonoid_core': Chem.MolFromSmarts('c1cc(c(cc1)O)-c1coc2c1ccc(=O)c(=O)c2'),  # Isoflavonoid skeleton
        'stilbene_core': Chem.MolFromSmarts('c1ccccc1C=Cc2ccccc2'),  # Stilbene scaffold
        'lignan_core': Chem.MolFromSmarts('c1cc(c(cc1)O)C[C@H](C2=CC=C(C=C2)O)O'),  # Lignan structure
    }

    # Check for matches to any of the patterns
    for name, pattern in patterns.items():
        if mol.HasSubstructMatch(pattern):
            return True, f"Contains {name.replace('_', ' ')} substructure"

    # Check for C6-C3 backbone (phenyl attached to 3-carbons)
    phenyl = Chem.MolFromSmarts('c1ccccc1')
    if mol.HasSubstructMatch(phenyl):
        # Find phenyl ring atoms
        phenyl_match = mol.GetSubstructMatch(phenyl)
        phenyl_atoms = set(phenyl_match)

        # Search for connected 3-carbon chain
        for atom_idx in phenyl_atoms:
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() not in phenyl_atoms:
                    # Check if neighbor is connected to a 3-carbon chain
                    chain = Chem.PathToSubmol(mol, [atom_idx, neighbor.GetIdx()])
                    if chain.GetNumHeavyAtoms() >= 3:
                        return True, "Contains phenylpropane skeleton"

    return False, "Does not contain phenylpropanoid substructure"

__metadata__ = {   
    'chemical_class': {   
        'name': 'phenylpropanoid',
        'definition': 'Any organic aromatic compound with a structure based on a phenylpropane skeleton. The class includes naturally occurring phenylpropanoid esters, flavonoids, anthocyanins, coumarins and many small phenolic molecules as well as their semi-synthetic and synthetic analogues. Phenylpropanoids are also precursors of lignin.',
        'parents': []},
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
        'test_proportion': 0.1},
    'message': None,
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'precision': None,
    'recall': None,
    'f1': None,
    'accuracy': None
}