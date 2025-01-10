"""
Classifies: CHEBI:33839 macromolecule
"""
"""
Classifies: macromolecule
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Descriptors
from rdkit.Chem import rdFMCS

def is_macromolecule(smiles: str):
    """
    Determines if a molecule is a macromolecule based on its SMILES string.
    A macromolecule is defined as a molecule of high relative molecular mass,
    comprising multiple repetitions of units derived from small molecules.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a macromolecule, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 1000:
        return False, f"Molecular weight is {mol_wt:.2f} Da, which is below the macromolecule threshold"

    # Attempt to detect repeating units
    # Create a list of fragments by breaking rotatable bonds
    frags = Chem.FragmentOnBonds(mol, mol.GetSubstructMatches(Chem.MolFromSmarts("[!$(*#*)&!D1]-!@[!$(*#*)&!D1]")), addDummies=False)
    frag_mols = Chem.GetMolFrags(frags, asMols=True)
    
    # Count unique fragments
    frag_smiles = [Chem.MolToSmiles(frag, isomericSmiles=True) for frag in frag_mols]
    unique_frags = set(frag_smiles)
    if len(unique_frags) < len(frag_smiles) / 2:
        return True, "Molecule has high molecular weight and contains repeating units"
    else:
        return False, "No significant repeating units detected despite high molecular weight"

    # If unable to determine, return None
    # return None, None

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'macromolecule',
        'definition': 'A macromolecule is a molecule of high relative molecular mass, the structure of which essentially comprises the multiple repetition of units derived, actually or conceptually, from molecules of low relative molecular mass.',
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
    'attempt': 0,
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