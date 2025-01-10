"""
Classifies: CHEBI:26658 sesquiterpenoid
"""
from rdkit import Chem

def is_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a sesquiterpenoid based on its SMILES string.
    Sesquiterpenoids are characterized by a C15 backbone that may be rearranged or modified,
    often retaining a sesquiterpene parent structure with additional functional groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a sesquiterpenoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbon atoms (flexible range to allow derivatives and complex forms)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 13 or c_count > 30:
        return False, f"Expected sesquiterpenoid-like carbon range, got {c_count}"

    # Define broad structural and functional patterns
    isoprene_unit = Chem.MolFromSmarts('C=C[C@H](C)C')  # basic isoprene pattern
    substructures = [
        Chem.MolFromSmarts('C1CCC(C)C1'),  # cycloalkane ring
        Chem.MolFromSmarts('[cH]1[cH][cH][cH][cH][cH]1'),  # aromatic-like ring
        Chem.MolFromSmarts('C=CO'),  # allylic alcohol
        Chem.MolFromSmarts('O=C(O)'),  # carboxylic acid or ester
        Chem.MolFromSmarts('C1=C[C@H](C)C[C@@H]1'),  # cycle with chirality
        Chem.MolFromSmarts('C=O'),  # carbonyl groups
        Chem.MolFromSmarts('COC'),  # methoxy groups
        Chem.MolFromSmarts('C#C'),  # acetylene group
    ]

    # Check for isoprene units
    isoprene_matches = mol.HasSubstructMatch(isoprene_unit)

    # Count matches for sesquiterpenoid characteristics
    matches = sum(mol.HasSubstructMatch(pattern) for pattern in substructures)

    # Empirically decide how many matches define a sesquiterpenoid
    if matches >= 3 and isoprene_matches:
        return True, "Contains multiple sesquiterpenoid characteristic features"

    return False, "Does not exhibit enough sesquiterpenoid characteristic features"