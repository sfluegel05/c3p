"""
Classifies: CHEBI:67142 nucleobase analogue
"""
"""
Classifies: nucleobase analogue
"""

from rdkit import Chem

def is_nucleobase_analogue(smiles: str):
    """
    Determines if a molecule is a nucleobase analogue based on its SMILES string.
    A nucleobase analogue is a molecule that can substitute for a normal nucleobase in nucleic acids.
    They typically have a purine or pyrimidine core with possible modifications,
    and can form hydrogen bonds similar to canonical nucleobases.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleobase analogue, False otherwise
        str: Reason for classification
    """
    from rdkit.Chem import AllChem

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Remove salts and keep only the largest fragment
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    mol = max(frags, default=mol, key=lambda m: m.GetNumAtoms())

    # Define SMARTS patterns for nucleobase cores (purine and pyrimidine analogues)
    nucleobase_patterns = [
        Chem.MolFromSmarts('c1cncnc1'),  # Pyrimidine core
        Chem.MolFromSmarts('n1cnc2c1ncnc2'),  # Purine core
        Chem.MolFromSmarts('c1ncnc1'),  # Pyrazine core
        Chem.MolFromSmarts('c1c[nH]cn1'),  # Imidazole-pyridine core
        Chem.MolFromSmarts('n1c(=O)[nH]c(=O)[nH]c1'),  # Uracil and analogues
        Chem.MolFromSmarts('O=c1cc[nH]c(=O)[nH]1'),  # Additional uracil-like pattern
        Chem.MolFromSmarts('n1c(=O)[nH]cc(=O)[nH]c1'),  # Cytosine and analogues
        Chem.MolFromSmarts('c1[nH]c2c(n1)ncnc2=O'),  # Guanine-like core
        Chem.MolFromSmarts('c1ncnc2ncnn12'),  # Aza-adenine core
        Chem.MolFromSmarts('c1ncnc2[nH]nnc12'),  # Triazolopyrimidine core
    ]

    # Define SMARTS patterns for sugar (ribose/deoxyribose) and phosphate groups
    sugar_patterns = [
        Chem.MolFromSmarts('C1OC([H])C([H])C1O'),  # Simplified sugar pattern
        Chem.MolFromSmarts('[C@H]1([O])[C@@H](O)[C@@H](O)[C@H](CO)O1'),  # Ribose
    ]
    phosphate_pattern = Chem.MolFromSmarts('O=P([O-])(O)O')  # Phosphate group

    # Exclude molecules with attached sugars or phosphates
    if any(mol.HasSubstructMatch(sugar) for sugar in sugar_patterns):
        return False, "Contains sugar group; likely a nucleoside or nucleotide"
    if mol.HasSubstructMatch(phosphate_pattern):
        return False, "Contains phosphate group; likely a nucleotide"

    # Check for nucleobase cores
    for pattern in nucleobase_patterns:
        if mol.HasSubstructMatch(pattern):
            # Exclude molecules with attached long chains or bulky groups
            num_atoms = mol.GetNumAtoms()
            num_rings = mol.GetRingInfo().NumRings()
            if num_atoms > 30:
                return False, "Molecule too large to be a simple nucleobase analogue"
            if num_rings > 3:
                return False, "Contains additional ring systems; not a nucleobase analogue"
            return True, "Contains nucleobase core structure"

    return False, "Does not contain nucleobase core structure"

__metadata__ = {  
    'chemical_class': {   
        'name': 'nucleobase analogue',
        'definition': 'A molecule that can substitute for a normal nucleobase in nucleic acids.',
    },
    'message': None,
    'success': True,
    'error': '',
}