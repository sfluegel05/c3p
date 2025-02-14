"""
Classifies: CHEBI:36976 nucleotide
"""
"""
Classifies: CHEBI:33561 nucleotide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_nucleotide(smiles: str):
    """
    Determines if a molecule is a nucleotide based on its SMILES string.
    A nucleotide is a nucleoside phosphate resulting from the condensation
    of the 3' or 5' hydroxy group of a nucleoside with phosphoric acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleotide, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns
    # Sugar (ribose or deoxyribose)
    sugar_pattern = Chem.MolFromSmarts("C1[C@H](O)[C@@H](O)[C@H](O)[C@H]1O")  # Furanose ring

    # Nitrogenous bases (purines and pyrimidines)
    purine_pattern = Chem.MolFromSmarts("c1ncnc2ncnn12")  # Purine base
    pyrimidine_pattern = Chem.MolFromSmarts("c1cncnc1=O")  # Pyrimidine base (uracil)
    cytosine_pattern = Chem.MolFromSmarts("c1cnc(N)nc1=O")  # Cytosine
    thymine_pattern = Chem.MolFromSmarts("Cc1cnc(=O)[nH]c1=O")  # Thymine
    guanine_pattern = Chem.MolFromSmarts("c1[nH]c2c(n1)nc(N)nc2=O")  # Guanine
    adenine_pattern = Chem.MolFromSmarts("c1[nH]c2c(n1)nc(N)nc2")  # Adenine

    # Phosphate group
    phosphate_pattern = Chem.MolFromSmarts("OP(=O)(O)[O-]")  # Phosphate group

    # Check for sugar moiety
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No ribose or deoxyribose sugar found"

    # Check for nitrogenous base
    base_found = False
    base_patterns = [purine_pattern, pyrimidine_pattern, cytosine_pattern,
                     thymine_pattern, guanine_pattern, adenine_pattern]
    for base_pattern in base_patterns:
        if mol.HasSubstructMatch(base_pattern):
            base_found = True
            break
    if not base_found:
        return False, "No nitrogenous base found"

    # Check for phosphate group
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"

    # Check for connection between sugar and base (glycosidic bond)
    # Define bond between sugar anomeric carbon and base nitrogen
    glycosidic_bond = Chem.MolFromSmarts("[C@H]1([O])[C@@H](O)[C@H](O)[C@H](O)1[*]")  # Sugar connected to base
    if not mol.HasSubstructMatch(glycosidic_bond):
        return False, "No glycosidic bond found between sugar and base"

    # Check for phosphate attached to 3' or 5' hydroxyl group of sugar
    # 5' phosphate pattern
    five_prime_phosphate = Chem.MolFromSmarts("O[P](=O)(O)OC[C@H]1O[C@@H](O)[C@H](O)[C@@H]1O")  # Phosphate on 5' carbon
    # 3' phosphate pattern
    three_prime_phosphate = Chem.MolFromSmarts("O[P](=O)(O)O[C@H]1O[C@@H](O)[C@H](O)[C@@H]1CO")  # Phosphate on 3' carbon

    if mol.HasSubstructMatch(five_prime_phosphate) or mol.HasSubstructMatch(three_prime_phosphate):
        return True, "Molecule is a nucleotide with phosphate group on sugar"
    else:
        return False, "Phosphate group not correctly attached to sugar"