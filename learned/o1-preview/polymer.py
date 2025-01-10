"""
Classifies: CHEBI:60027 polymer
"""
"""
Classifies: polymer
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polymer(smiles: str):
    """
    Determines if a molecule is a polymer based on its SMILES string.
    A polymer is characterized by high molecular weight and large size due to repeating units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polymer, False otherwise
        str: Reason for classification
    """

    # Parse SMILES with error handling
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=True)
        if mol is None:
            return False, "Invalid SMILES string"
    except Chem.rdchem.KekulizeException as e:
        return False, f"SMILES Parsing Error: {str(e)}"
    except Exception as e:
        return False, f"SMILES Parsing Error: {str(e)}"

    # Remove counter ions and small fragments, keep the largest fragment
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if not frags:
        return False, "No valid molecular fragments found"
    elif len(frags) > 1:
        mol = max(frags, key=lambda m: m.GetNumAtoms())
    else:
        mol = frags[0]

    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    # Count number of atoms
    num_atoms = mol.GetNumAtoms()
    # Count number of rotatable bonds
    num_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)

    # Set threshold values (these can be adjusted based on domain knowledge)
    MW_THRESHOLD = 1000  # Molecular weight threshold in Daltons
    ATOM_THRESHOLD = 100 # Number of atoms threshold
    ROTATABLE_BOND_THRESHOLD = 20  # Number of rotatable bonds

    if mol_wt > MW_THRESHOLD and num_atoms > ATOM_THRESHOLD and num_rotatable > ROTATABLE_BOND_THRESHOLD:
        return True, f"Molecule has high molecular weight ({mol_wt:.2f} Da), {num_atoms} atoms, and {num_rotatable} rotatable bonds, indicative of a polymer"
    else:
        reasons = []
        if mol_wt <= MW_THRESHOLD:
            reasons.append(f"Molecular weight ({mol_wt:.2f} Da) is below threshold ({MW_THRESHOLD} Da)")
        if num_atoms <= ATOM_THRESHOLD:
            reasons.append(f"Number of atoms ({num_atoms}) is below threshold ({ATOM_THRESHOLD})")
        if num_rotatable <= ROTATABLE_BOND_THRESHOLD:
            reasons.append(f"Number of rotatable bonds ({num_rotatable}) is below threshold ({ROTATABLE_BOND_THRESHOLD})")
        reason_str = "; ".join(reasons)
        return False, reason_str


__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'polymer',
        'definition': 'A polymer is a mixture, which is composed of macromolecules of different kinds and which may be differentiated by composition, length, degree of branching etc.',
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
    'attempt': 1,
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