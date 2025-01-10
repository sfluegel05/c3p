"""
Classifies: CHEBI:46662 mineral
"""
"""
Classifies: CHEBI:46662 mineral
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_mineral(smiles: str):
    """
    Determines if a molecule is a mineral based on its SMILES string.
    A mineral is a naturally occurring inorganic substance formed through geological processes.
    This function checks for the absence of carbon-hydrogen bonds and common organic functional groups,
    and the presence of metal ions combined with inorganic anions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a mineral, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Exclude molecules with aromatic rings (typical of organic compounds)
    if mol.GetAromaticAtomCount() > 0:
        return False, "Contains aromatic rings typical of organic molecules"

    # Exclude molecules with carbon-hydrogen bonds (indicative of organic compounds)
    for bond in mol.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        if (atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 1) or \
           (atom1.GetAtomicNum() == 1 and atom2.GetAtomicNum() == 6):
            return False, "Contains carbon-hydrogen bonds typical of organic molecules"

    # Exclude molecules with typical organic functional groups
    organic_functional_groups = [
        Chem.MolFromSmarts("[CX3]=[OX1]"),                 # Carbonyl group
        Chem.MolFromSmarts("[CX3]=[CX3]"),                 # Alkene
        Chem.MolFromSmarts("[CX2]#[CX1]"),                 # Alkyne
        Chem.MolFromSmarts("[#6][OX2H]"),                  # Alcohol
        Chem.MolFromSmarts("[#6][NX3;H2,H1;!$([N][O])]"),  # Amine
        Chem.MolFromSmarts("[#6][#16][#6]"),               # Thioether
        Chem.MolFromSmarts("[CX3](=O)[OX2H1]"),            # Carboxylic acid
        Chem.MolFromSmarts("c"),                           # Aromatic carbons
    ]
    for fg in organic_functional_groups:
        if mol.HasSubstructMatch(fg):
            return False, "Contains functional groups typical of organic molecules"

    # Identify metals in the molecule
    metals = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() in rdchem.PeriodicTable.GetElementSymbolList()
              and atom.GetAtomicNum() != 6 and atom.GetAtomicNum() != 1]
    if not metals:
        return False, "No metal atoms found in the molecule"

    # Check for inorganic anions
    inorganic_anions = [
        Chem.MolFromSmarts("[O-2]"),            # Oxide
        Chem.MolFromSmarts("[S-2]"),            # Sulfide
        Chem.MolFromSmarts("[N-3]"),            # Nitride
        Chem.MolFromSmarts("[F-]"),             # Fluoride
        Chem.MolFromSmarts("[Cl-]"),            # Chloride
        Chem.MolFromSmarts("[Br-]"),            # Bromide
        Chem.MolFromSmarts("[I-]"),             # Iodide
        Chem.MolFromSmarts("[OH-]"),            # Hydroxide
        Chem.MolFromSmarts("O=[SX4](=O)[O-]"),  # Sulfate
        Chem.MolFromSmarts("O=[NX3](=O)[O-]"),  # Nitrate
        Chem.MolFromSmarts("[CO3]"),            # Carbonate
        Chem.MolFromSmarts("P(=O)([O-])([O-])[O-]"),  # Phosphate
        Chem.MolFromSmarts("[Si]([O-])([O-])[O-]"),   # Silicate
    ]
    anion_found = False
    for anion in inorganic_anions:
        if mol.HasSubstructMatch(anion):
            anion_found = True
            break
    if not anion_found:
        return False, "No inorganic anions found in the molecule"

    # Exclude organometallics (metal-carbon bonds)
    for bond in mol.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        if (atom1.GetAtomicNum() in range(21, 31) + range(39, 49) + range(72, 81) and atom2.GetAtomicNum() == 6) or \
           (atom2.GetAtomicNum() in range(21, 31) + range(39, 49) + range(72, 81) and atom1.GetAtomicNum() == 6):
            return False, "Contains metal-carbon bond typical of organometallic compounds"

    return True, "Molecule meets criteria for mineral classification"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:46662',
        'name': 'mineral',
        'definition': 'In general, a mineral is a chemical substance that is normally crystalline formed and has been formed as a result of geological processes. The term also includes metamict substances (naturally occurring, formerly crystalline substances whose crystallinity has been destroyed by ionising radiation) and can include naturally occurring amorphous substances that have never been crystalline (\'mineraloids\') such as georgite and calciouranoite as well as substances formed by the action of geological processes on biogenic compounds (\'biogenic minerals\').',
        'parents': ['CHEBI:24437']
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
    'attempt': 3,
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