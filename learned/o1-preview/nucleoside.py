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
    A nucleoside consists of a nucleobase attached to a ribose or deoxyribose sugar via an N-glycosidic bond.

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

    # Define SMARTS patterns for nucleobases (more general)
    nucleobase_smarts = [
        # Purines (adenine, guanine, xanthine, hypoxanthine)
        '[nH]1c2[nH]cnc2nc1',           # Adenine-like
        'n1(cnc2c1ncn2)[NH2]',          # Guanine-like
        'O=C1NC(=O)c2[nH]cnc12',        # Xanthine-like
        'O=C1NC=NC2=C1N=CN2',           # Hypoxanthine-like
        # Pyrimidines (cytosine, thymine, uracil)
        'C1=NC(=O)NC(=O)C1',            # Uracil-like
        'C1=NC(=O)NC(=O)C1C',           # Thymine-like
        'C1=C(N)NC=NC1=O',              # Cytosine-like
    ]

    nucleobase_patterns = [Chem.MolFromSmarts(smarts) for smarts in nucleobase_smarts]

    # Define SMARTS patterns for ribose and deoxyribose sugars (more general)
    ribose_smarts = '[C@H]1([O])[C@@H]([O])[C@H]([O])[C@@H](CO)[O1]'  # Ribose ring
    deoxyribose_smarts = '[C@H]1([O])[C@@H]([O])[C@H](CO)[C@@H](CO)[O1]'  # Deoxyribose ring

    ribose_pattern = Chem.MolFromSmarts(ribose_smarts)
    deoxyribose_pattern = Chem.MolFromSmarts(deoxyribose_smarts)

    # Check for nucleobase presence
    has_nucleobase = False
    for pattern in nucleobase_patterns:
        if mol.HasSubstructMatch(pattern):
            has_nucleobase = True
            break
    if not has_nucleobase:
        return False, "No nucleobase found"

    # Check for sugar moiety presence
    has_ribose = mol.HasSubstructMatch(ribose_pattern)
    has_deoxyribose = mol.HasSubstructMatch(deoxyribose_pattern)
    if not (has_ribose or has_deoxyribose):
        return False, "No ribose or deoxyribose sugar moiety found"

    # Check for N-glycosidic bond between sugar and nucleobase
    # Find nitrogen atoms in nucleobase
    nucleobase_nitrogens = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7 and any(
        mol.HasSubstructMatch(Chem.MolFromSmarts(f"[#{atom.GetAtomicNum()}]"), useChirality=False))]

    # Find anomeric carbon in sugar (the carbon connected to the nucleobase)
    if has_ribose:
        sugar_pattern = ribose_pattern
    else:
        sugar_pattern = deoxyribose_pattern

    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if not sugar_matches:
        return False, "Sugar moiety not properly matched"

    # Assume first match is the sugar
    sugar_atoms = sugar_matches[0]
    sugar_atom_set = set(sugar_atoms)

    # Find the anomeric carbon (CA) in the sugar ring
    anomeric_carbons = [idx for idx in sugar_atoms if mol.GetAtomWithIdx(idx).GetDegree() == 3]
    if not anomeric_carbons:
        # Fallback: Find carbon in sugar connected to oxygen outside the ring
        for idx in sugar_atoms:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 6:
                for nbr in atom.GetNeighbors():
                    if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in sugar_atom_set:
                        anomeric_carbons.append(idx)
                        break
    if not anomeric_carbons:
        return False, "Anomeric carbon not found in sugar"

    # Check for bond between anomeric carbon and nucleobase nitrogen
    glycosidic_bond_found = False
    for ac_idx in anomeric_carbons:
        anomeric_carbon = mol.GetAtomWithIdx(ac_idx)
        for nbr in anomeric_carbon.GetNeighbors():
            if nbr.GetIdx() in nucleobase_nitrogens:
                glycosidic_bond_found = True
                break
        if glycosidic_bond_found:
            break

    if not glycosidic_bond_found:
        return False, "No N-glycosidic bond connecting nucleobase and sugar found"

    return True, "Molecule is a nucleoside with nucleobase attached to sugar via N-glycosidic bond"