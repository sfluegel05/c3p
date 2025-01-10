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
    They typically have a core structure similar to canonical nucleobases (adenine, guanine, cytosine, thymine, uracil)
    with possible modifications, but without attached sugar or phosphate groups.

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

    # Define SMARTS patterns for nucleobase analogues
    nucleobase_patterns = {
        'adenine_like': Chem.MolFromSmarts('n1c([#6,#7,#8,#9,#16])[nH]c2c1ncnc2'),
        'guanine_like': Chem.MolFromSmarts('n1c([#6,#7,#8,#9,#16])nc2c1[nH]c(=O)[nH]c2'),
        'cytosine_like': Chem.MolFromSmarts('n1c([#6,#7,#8,#9,#16])nc([#6,#7,#8,#9,#16])cc1=O'),
        'uracil_like': Chem.MolFromSmarts('O=C1NC(=O)C([#6,#7,#8,#9,#16])=CN1'),
        'thymine_like': Chem.MolFromSmarts('O=C1NC(=O)C([#6,#7,#8,#9,#16])=C(N1)[#6]')
    }

    # Define SMARTS patterns for sugar (ribose/deoxyribose)
    sugar_patterns = [
        Chem.MolFromSmarts('C1OC([H])C([H])C1O'),  # Simplified sugar pattern
        Chem.MolFromSmarts('[C@H]1([O])[C@@H](O)[C@@H](O)[C@H](CO)O1'),  # Ribose
        Chem.MolFromSmarts('[C@H]1([O])[C@@H](O)[C@@H](O)[C@H](CO)O1')   # Deoxyribose
    ]

    # Define SMARTS pattern for phosphate group
    phosphate_pattern = Chem.MolFromSmarts('P(=O)(O)O')

    # Check for sugar or phosphate groups
    for sugar_pattern in sugar_patterns:
        if mol.HasSubstructMatch(sugar_pattern):
            return False, "Contains sugar group; likely a nucleoside or nucleotide"
    if mol.HasSubstructMatch(phosphate_pattern):
        return False, "Contains phosphate group; likely a nucleotide"

    # Check for nucleobase core structures allowing modifications
    for name, pattern in nucleobase_patterns.items():
        if mol.HasSubstructMatch(pattern):
            return True, f"Contains nucleobase analogue core structure similar to {name}"

    # Exclude molecules with additional large ring systems
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings > 2:
        return False, "Contains additional ring systems; not a nucleobase analogue"

    # Check molecular weight - nucleobases are usually less than 200 Da
    mol_wt = Chem.rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 300:
        return False, "Molecular weight too high for a nucleobase analogue"

    # Check for attached long chains (e.g., peptides or large substituents)
    if any(atom.GetDegree() > 3 for atom in mol.GetAtoms()):
        return False, "Contains bulky substituents; unlikely to be a nucleobase analogue"

    # As a fallback, check for pyrimidine or purine core without sugar/phosphate
    pyrimidine_core = Chem.MolFromSmarts('c1cncnc1')
    purine_core = Chem.MolFromSmarts('n1cnc2c1ncnc2')
    if mol.HasSubstructMatch(pyrimidine_core) or mol.HasSubstructMatch(purine_core):
        return True, "Contains purine or pyrimidine core; possible nucleobase analogue"

    return False, "Does not contain nucleobase analogue core structure"