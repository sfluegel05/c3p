"""
Classifies: CHEBI:33838 nucleoside
"""
"""
Classifies: CHEBI:33838 nucleoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_nucleoside(smiles: str):
    """
    Determines if a molecule is a nucleoside based on its SMILES string.
    A nucleoside is an N-glycosyl compound that has both a nucleobase
    (normally adenine, guanine, xanthine, thymine, cytosine, or uracil)
    and either a ribose or deoxyribose sugar as functional parents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleoside, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Remove salts and counterions
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if len(frags) > 1:
        mol = max(frags, default=mol, key=lambda m: m.GetNumAtoms())

    # Check for phosphate groups to exclude nucleotides
    phosphate_smarts = Chem.MolFromSmarts('P(=O)(O)(O)O')
    if mol.HasSubstructMatch(phosphate_smarts):
        return False, "Contains phosphate group, likely a nucleotide"

    # Define SMARTS patterns for ribose and deoxyribose
    ribose_smarts = Chem.MolFromSmarts('O[C@H]1[C@@H](O)[C@H](O)[C@@H](CO)O1')
    deoxyribose_smarts = Chem.MolFromSmarts('O[C@H]1[C@@H](O)[C@H](CO)[C@@H]([CH2])O1')

    # Check for sugar moiety
    has_sugar = mol.HasSubstructMatch(ribose_smarts) or mol.HasSubstructMatch(deoxyribose_smarts)
    if not has_sugar:
        return False, "No ribose or deoxyribose sugar found"

    # Define SMARTS patterns for common nucleobases
    nucleobase_smarts = [
        # Adenine
        'n1c[nH]c2c1ncn2',
        # Guanine
        'n1c2c(c(=O)[nH]1)nc[nH]2',
        # Cytosine
        'n1c(N)ccn1',
        # Uracil/Thymine
        'O=C1NC(=O)C=C1',
        # Xanthine
        'O=C1NC(=O)Nc2ncnc12',
    ]

    # Check for nucleobase
    has_base = False
    for base in nucleobase_smarts:
        base_pattern = Chem.MolFromSmarts(base)
        if mol.HasSubstructMatch(base_pattern):
            has_base = True
            break

    if not has_base:
        return False, "No common nucleobase found"

    # Check for N-glycosidic bond between sugar and base
    glycosidic_bond_smarts = Chem.MolFromSmarts('[C@H]1([O])[C@@H]([O])[C@H]([O])[C@@H](CO1)N2C=NC=NC2')  # Generic pattern
    if not mol.HasSubstructMatch(glycosidic_bond_smarts):
        return False, "No N-glycosidic bond between sugar and nucleobase found"

    return True, "Contains nucleobase attached to sugar via N-glycosidic bond"