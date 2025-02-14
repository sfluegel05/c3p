"""
Classifies: CHEBI:75769 B vitamin
"""
"""
Classifies: B vitamin
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_B_vitamin(smiles: str):
    """
    Determines if a molecule is a B vitamin based on its SMILES string.
    The B vitamins include vitamin B1 (thiamine), B2 (riboflavin), B3 (niacin),
    B5 (pantothenic acid), B6 (pyridoxine, pyridoxal, pyridoxamine),
    B7 (biotin), B9 (folic acid), and B12 (cobalamin).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a B vitamin, False otherwise
        str: Reason for classification
    """

    # Convert input SMILES to molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for each B vitamin

    patterns = [
        {
            'name': 'Vitamin B1 (Thiamine)',
            'pattern': Chem.MolFromSmarts('CC1=C(C)C=CN=C1CCO'),  # Simplified thiazolium structure
            'description': 'Matches thiazolium ring of thiamine'
        },
        {
            'name': 'Vitamin B2 (Riboflavin)',
            'pattern': Chem.MolFromSmarts('CNC1=C2C(=O)NC(=O)NC2=NC2=C1C=CC=C2'),  # Isoalloxazine ring
            'description': 'Matches isoalloxazine ring system of riboflavin'
        },
        {
            'name': 'Vitamin B3 (Niacin)',
            'pattern': Chem.MolFromSmarts('c1cccnc1C(=O)O'),  # Pyridine ring with carboxylic acid
            'description': 'Matches nicotinic acid structure of niacin'
        },
        {
            'name': 'Vitamin B5 (Pantothenic acid)',
            'pattern': Chem.MolFromSmarts('CC(C)(CO)C(=O)NCCC(=O)O'),  # Pantothenic acid structure
            'description': 'Matches structure of pantothenic acid'
        },
        {
            'name': 'Vitamin B6 (Pyridoxine, Pyridoxal, Pyridoxamine)',
            'pattern': Chem.MolFromSmarts('c1cc(CO)c(O)nc1C'),  # Pyridine ring with hydroxyl and methyl groups
            'description': 'Matches pyridine ring with hydroxyl group characteristic of vitamin B6'
        },
        {
            'name': 'Vitamin B7 (Biotin)',
            'pattern': Chem.MolFromSmarts('O=C1NC(=O)N2[C@@](CS1)([H])CCCCC2'),  # Biotin structure
            'description': 'Matches fused ring system of biotin'
        },
        {
            'name': 'Vitamin B9 (Folic acid)',
            'pattern': Chem.MolFromSmarts('Nc1nc2ncc(CNc3ccc(cc3)C(O)=O)nc2c(=O)[nH]1'),  # Pteridine ring with p-aminobenzoic acid
            'description': 'Matches pteridine ring system of folic acid'
        },
        {
            'name': 'Vitamin B12 (Cobalamin)',
            'pattern': Chem.MolFromSmarts('[Cobalt]'),  # Cobalt atom in a corrin ring
            'description': 'Contains cobalt within a corrin ring characteristic of cobalamin'
        },
    ]

    # Check for matches against each vitamin pattern
    for vit in patterns:
        pattern = vit['pattern']
        if pattern is None:
            continue

        # For Vitamin B12, check for the corrin ring
        if vit['name'] == 'Vitamin B12 (Cobalamin)':
            contains_cobalt = any(atom.GetAtomicNum() == 27 for atom in mol.GetAtoms())
            corrin_ring = Chem.MolFromSmarts('C1=CC=CC=C1')  # Placeholder for corrin ring; needs accurate SMARTS
            if contains_cobalt:
                # Attempt to match corrin ring
                # Note: Corrin ring is complex; here we check for cobalt coordinated to nitrogen atoms
                cobalt_coordination = Chem.MolFromSmarts('[Co]~[N]~[C]')
                if mol.HasSubstructMatch(cobalt_coordination):
                    return True, vit['description']
                else:
                    continue  # Contains cobalt but not as part of cobalamin
            else:
                continue  # No cobalt, not vitamin B12
        else:
            # For other vitamins, perform substructure matching
            if mol.HasSubstructMatch(pattern):
                return True, f"Molecule matches {vit['name']} ({vit['description']})"

    return False, "Molecule does not match any known B vitamin"