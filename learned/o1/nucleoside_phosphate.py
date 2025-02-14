"""
Classifies: CHEBI:25608 nucleoside phosphate
"""
"""
Classifies: nucleoside phosphate
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nucleoside_phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside phosphate based on its SMILES string.
    A nucleoside phosphate is a nucleobase-containing molecule where one or more of
    the sugar hydroxy groups has been converted into a mono- or poly-phosphate.
    This includes both nucleotides and non-nucleotide nucleoside phosphates.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is a nucleoside phosphate, False otherwise.
        str: Reason for classification.
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define nucleobase patterns (common nucleobases and modifications)
    nucleobase_smarts_list = [
        # Purine base (adenine, guanine, hypoxanthine)
        'n1cnc2ncnc(N)c12',    # Adenine
        'NC1=NC2=NC=NC(N)=C2N1',  # Guanine
        'n1cnc2nc[nH]c2n1',    # Hypoxanthine
        # Pyrimidine base (cytosine, uracil, thymine)
        'n1cc(N)cnc1=O',       # Cytosine
        'O=C1C=CN(C)C=C1',     # Thymine
        'O=C1C=CC(N)=CN1',     # Uracil
        # Modified bases
        'n1c(=O)[nH]c2c1ncnc2',  # Xanthine
        'n1c(=O)c[nH]c1=O',    # Alloxazine (in flavins)
        # General purine and pyrimidine patterns
        'c1ncnc2c1ncn2',       # Purine skeleton
        'c1[nH]cnc1',          # Pyrimidine skeleton
    ]
    nucleobase_patterns = [Chem.MolFromSmarts(smarts) for smarts in nucleobase_smarts_list]

    # Search for nucleobase substructure
    nucleobase_found = False
    for pattern in nucleobase_patterns:
        if mol.HasSubstructMatch(pattern):
            nucleobase_found = True
            break
    if not nucleobase_found:
        return False, "No nucleobase found"

    # Define flexible sugar pattern (ribose or deoxyribose with possible modifications)
    sugar_smarts = '[C@H1]([O])[C@@H]([O])[C@H]([O])[C@@H]([O])[CH2]'
    sugar_pattern = Chem.MolFromSmarts(sugar_smarts)

    # Search for sugar substructure
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if not sugar_matches:
        return False, "No sugar moiety found"

    # Check for glycosidic bond between nucleobase and sugar
    glycosidic_bond_found = False
    for sugar_match in sugar_matches:
        sugar_atoms = set(sugar_match)
        for atom_idx in sugar_atoms:
            atom = mol.GetAtomWithIdx(atom_idx)
            for bond in atom.GetBonds():
                neighbor = bond.GetOtherAtom(atom)
                # Check if neighbor is part of the nucleobase
                for pattern in nucleobase_patterns:
                    if mol.GetSubstructMatch(pattern):
                        if neighbor.GetIdx() in mol.GetSubstructMatch(pattern):
                            glycosidic_bond_found = True
                            break
                if glycosidic_bond_found:
                    break
            if glycosidic_bond_found:
                break
        if glycosidic_bond_found:
            break
    if not glycosidic_bond_found:
        return False, "No glycosidic bond between nucleobase and sugar found"

    # Define phosphate group pattern (mono- or polyphosphate)
    phosphate_smarts = 'OP(=O)(O)O'  # Phosphate group
    phosphate_pattern = Chem.MolFromSmarts(phosphate_smarts)

    # Search for phosphate groups
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if not phosphate_matches:
        return False, "No phosphate group found"

    # Check if phosphate is attached to sugar hydroxyl groups
    phosphate_attached_to_sugar = False
    for phosphate_match in phosphate_matches:
        phosphate_atom_indices = set(phosphate_match)
        for sugar_match in sugar_matches:
            sugar_atom_indices = set(sugar_match)
            for atom_idx in phosphate_atom_indices:
                atom = mol.GetAtomWithIdx(atom_idx)
                for bond in atom.GetBonds():
                    neighbor = bond.GetOtherAtom(atom)
                    if neighbor.GetIdx() in sugar_atom_indices:
                        phosphate_attached_to_sugar = True
                        break
                if phosphate_attached_to_sugar:
                    break
            if phosphate_attached_to_sugar:
                break
        if phosphate_attached_to_sugar:
            break
    if not phosphate_attached_to_sugar:
        return False, "No phosphate group attached to sugar"

    return True, "Molecule is a nucleoside phosphate with nucleobase, sugar, and phosphate group(s)"