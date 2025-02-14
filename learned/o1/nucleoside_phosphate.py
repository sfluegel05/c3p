"""
Classifies: CHEBI:25608 nucleoside phosphate
"""
"""
Classifies: nucleoside phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_nucleoside_phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside phosphate based on its SMILES string.
    A nucleoside phosphate is a nucleobase-containing molecule where one or more of
    the sugar hydroxy groups has been converted into a mono- or poly-phosphate.

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

    # Define nucleobase patterns
    adenine_smarts = "n1cnc2ncnc12"
    guanine_smarts = "n1c(=O)[nH]c2ncnc2c1=O"
    cytosine_smarts = "n1c(=O)nc[nH]c1"
    thymine_smarts = "Cn1c(=O)[nH]c(=O)c1"
    uracil_smarts = "O=C1NC=CC(=O)N1"
    nucleobase_patterns = [
        Chem.MolFromSmarts(adenine_smarts),
        Chem.MolFromSmarts(guanine_smarts),
        Chem.MolFromSmarts(cytosine_smarts),
        Chem.MolFromSmarts(thymine_smarts),
        Chem.MolFromSmarts(uracil_smarts)
    ]

    # Check for presence of nucleobase
    nucleobase_found = False
    for base_pattern in nucleobase_patterns:
        if mol.HasSubstructMatch(base_pattern):
            nucleobase_found = True
            break
    if not nucleobase_found:
        return False, "No nucleobase found"

    # Define sugar moiety pattern (ribose or deoxyribose)
    sugar_pattern = Chem.MolFromSmarts("C1OC[C@H](O)[C@@H]1O")  # Simplified ribose ring
    deoxysugar_pattern = Chem.MolFromSmarts("C1OC[C@H](O)[C@@H]1")  # Deoxyribose lacks one OH

    # Check for sugar moiety
    if not mol.HasSubstructMatch(sugar_pattern) and not mol.HasSubstructMatch(deoxysugar_pattern):
        return False, "No ribose or deoxyribose sugar moiety found"

    # Check for glycosidic bond between nucleobase and sugar
    glycosidic_bond_found = False
    # Get the index of atoms matching the sugar pattern
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    deoxysugar_matches = mol.GetSubstructMatches(deoxysugar_pattern)
    sugar_matches += deoxysugar_matches  # Combine matches
    nucleobase_matches = []
    for base_pattern in nucleobase_patterns:
        matches = mol.GetSubstructMatches(base_pattern)
        if matches:
            nucleobase_matches.extend(matches)

    if not nucleobase_matches:
        return False, "Nucleobase pattern not found in molecule"

    # Check for glycosidic bond between nucleobase and sugar
    for sugar_match in sugar_matches:
        for nucleobase_match in nucleobase_matches:
            for sugar_atom_idx in sugar_match:
                sugar_atom = mol.GetAtomWithIdx(sugar_atom_idx)
                if sugar_atom.GetAtomicNum() == 6:  # Carbon
                    # Check if this carbon is connected to nucleobase nitrogen
                    for bond in sugar_atom.GetBonds():
                        neighbor = bond.GetOtherAtom(sugar_atom)
                        neighbor_idx = neighbor.GetIdx()
                        if neighbor_idx in nucleobase_match and neighbor.GetAtomicNum() == 7:
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

    # Define phosphate group pattern attached to sugar's hydroxyl group
    phosphate_pattern = Chem.MolFromSmarts("OP(=O)(O)O")  # Phosphate group

    # Check for phosphate groups attached to sugar
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if not phosphate_matches:
        return False, "No phosphate groups attached to the sugar"

    phosphate_found = False
    for phosphate_match in phosphate_matches:
        phosphorus_idx = phosphate_match[1]  # Index of P atom
        phosphorus_atom = mol.GetAtomWithIdx(phosphorus_idx)
        # Check if phosphate is attached to sugar oxygen
        for bond in phosphorus_atom.GetBonds():
            neighbor = bond.GetOtherAtom(phosphorus_atom)
            if neighbor.GetAtomicNum() == 8:  # Oxygen
                for nb_bond in neighbor.GetBonds():
                    nb_neighbor = nb_bond.GetOtherAtom(neighbor)
                    if nb_neighbor.GetAtomicNum() == 6 and nb_neighbor.IsInRing():  # Carbon in ring (sugar)
                        phosphate_found = True
                        break
                if phosphate_found:
                    break
            if phosphate_found:
                break
        if phosphate_found:
            break
    if not phosphate_found:
        return False, "No phosphate group attached to sugar"

    return True, "Molecule is a nucleoside phosphate with nucleobase, sugar, and phosphate group(s)"