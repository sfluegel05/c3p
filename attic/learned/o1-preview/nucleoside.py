"""
Classifies: CHEBI:33838 nucleoside
"""
"""
Classifies: nucleoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_nucleoside(smiles: str):
    """
    Determines if a molecule is a nucleoside based on its SMILES string.
    A nucleoside is composed of a nucleobase attached to a ribose or deoxyribose sugar via an N-glycosidic bond.

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

    # Define SMARTS patterns for nucleobases
    nucleobase_smarts = {
        'adenine': 'n1cnc2c1ncnc2',
        'guanine': 'n1c(nc2c1ncnc2O)N',
        'cytosine': 'N1C=CC(=O)NC1=O',
        'thymine': 'C1=CN(C(=O)NC1=O)C',
        'uracil': 'O=C1NC(=O)C=C1',
        'xanthine': 'O=C1NC(=O)c2[nH]cnc12',
    }

    nucleobase_patterns = {name: Chem.MolFromSmarts(smarts) for name, smarts in nucleobase_smarts.items()}

    # Define SMARTS patterns for ribose and deoxyribose sugars
    ribose_smarts = '[C@H]1(O)[C@H](O)[C@@H](O)[C@H](O)[C@@H]1O'  # Ribose
    deoxyribose_smarts = '[C@H]1(O)[C@H](O)[C@@H](O)[C@H](O)[C@@H]1O'  # Deoxyribose (adjusted to allow for deoxy positions)

    ribose_pattern = Chem.MolFromSmarts(ribose_smarts)
    deoxyribose_pattern = Chem.MolFromSmarts(deoxyribose_smarts)

    # Check if nucleobase is present
    has_nucleobase = False
    for name, pattern in nucleobase_patterns.items():
        if mol.HasSubstructMatch(pattern):
            has_nucleobase = True
            base_name = name
            break
    if not has_nucleobase:
        return False, "No nucleobase found"

    # Check if sugar moiety is present
    has_ribose = mol.HasSubstructMatch(ribose_pattern)
    has_deoxyribose = mol.HasSubstructMatch(deoxyribose_pattern)

    if not (has_ribose or has_deoxyribose):
        return False, "No ribose or deoxyribose sugar moiety found"

    # Identify the sugar pattern used
    sugar_pattern = ribose_pattern if has_ribose else deoxyribose_pattern

    # Find the anomeric carbon in the sugar (carbon connected to two oxygen atoms)
    sugar_match = mol.GetSubstructMatch(sugar_pattern)
    if not sugar_match:
        return False, "Sugar moiety not properly matched"

    sugar_atoms = [mol.GetAtomWithIdx(idx) for idx in sugar_match]
    anomeric_carbons = [atom for atom in sugar_atoms if atom.GetAtomicNum() == 6 and
                        sum(1 for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 8) == 2]
    if not anomeric_carbons:
        return False, "Anomeric carbon not found in sugar"

    anomeric_carbon = anomeric_carbons[0]

    # Find nitrogen atoms in the nucleobase that could be linked to the sugar
    nucleobase_match = mol.GetSubstructMatch(nucleobase_patterns[base_name])
    nucleobase_atoms = [mol.GetAtomWithIdx(idx) for idx in nucleobase_match]
    nitrogen_atoms = [atom for atom in nucleobase_atoms if atom.GetAtomicNum() == 7]

    # Check for N-glycosidic bond between anomeric carbon and nucleobase nitrogen
    glycosidic_bond_found = False
    for bond in anomeric_carbon.GetBonds():
        neighbor = bond.GetOtherAtom(anomeric_carbon)
        if neighbor in nitrogen_atoms:
            glycosidic_bond_found = True
            break

    if not glycosidic_bond_found:
        return False, "No N-glycosidic bond connecting nucleobase and sugar found"

    return True, f"Contains {base_name} nucleobase attached to {'ribose' if has_ribose else 'deoxyribose'} sugar via N-glycosidic bond"