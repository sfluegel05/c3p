"""
Classifies: CHEBI:51129 nitrohydrocarbon
"""
"""
Classifies: CHEBI:XXXX nitrohydrocarbon
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_nitrohydrocarbon(smiles: str):
    """
    Determines if a molecule is a nitrohydrocarbon based on its SMILES string.
    A nitrohydrocarbon is a hydrocarbon in which one or more hydrogens has been replaced by nitro groups (C–NO2).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nitrohydrocarbon, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for nitro groups attached to carbon
    nitro_pattern = Chem.MolFromSmarts('[CX3](N(=O)=O)')  # Carbon attached to nitro group
    nitro_matches = mol.GetSubstructMatches(nitro_pattern)
    if not nitro_matches:
        return False, "No nitro groups attached to carbon found"

    # Clone molecule to remove nitro groups
    mol_no_nitro = Chem.RWMol(mol)
    atoms_to_remove = []
    for match in nitro_matches:
        carbon_idx = match[0]
        nitrogen_idx = match[1]
        oxygen1_idx = match[2]
        oxygen2_idx = match[3]
        atoms_to_remove.extend([nitrogen_idx, oxygen1_idx, oxygen2_idx])
        # Change the bond between carbon and nitrogen to single bond and to hydrogen
        mol_no_nitro.GetAtomWithIdx(carbon_idx).SetNumExplicitHs(mol_no_nitro.GetAtomWithIdx(carbon_idx).GetTotalNumHs() + 1)
        mol_no_nitro.GetAtomWithIdx(carbon_idx).UpdatePropertyCache()

    # Remove atoms in reverse order to maintain indices
    for idx in sorted(set(atoms_to_remove), reverse=True):
        mol_no_nitro.RemoveAtom(idx)

    # Update molecule
    mol_no_nitro.UpdatePropertyCache(strict=False)

    # Check that the remaining molecule contains only carbon and hydrogen atoms
    for atom in mol_no_nitro.GetAtoms():
        if atom.GetAtomicNum() not in (1, 6):  # H, C
            return False, f"Hydrocarbon part contains atom other than C and H: {atom.GetSymbol()}"

    # Check that there are no other heteroatoms in the original molecule besides nitro groups
    allowed_atomic_nums = {1, 6, 7, 8}  # H, C, N, O
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, f"Molecule contains heteroatom other than nitro group: {atom.GetSymbol()}"

    return True, "Molecule is a nitrohydrocarbon with nitro groups attached to a hydrocarbon"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:XXXX',
        'name': 'nitrohydrocarbon',
        'definition': 'A C-nitro compound that is a hydrocarbon in which one or more of the hydrogens has been replaced by nitro groups.',
        'parents': ['CHEBI:35474', 'CHEBI:62834']
    },
    'config': {
        'llm_model_name': 'INSERT_MODEL_NAME_HERE',
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