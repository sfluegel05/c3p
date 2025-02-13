"""
Classifies: CHEBI:36976 nucleotide
"""
"""
Classifies: CHEBI:60510 nucleotide
A nucleotide is a nucleoside phosphate resulting from the condensation of the 3 or 5 hydroxy group of a nucleoside with phosphoric acid.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_nucleotide(smiles: str):
    """
    Determines if a molecule is a nucleotide based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleotide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define nucleobase SMARTS patterns
    adenine = Chem.MolFromSmarts("nc1ncnc2n(cnc12)[CH]")
    guanine = Chem.MolFromSmarts("nc1nc2n(cnc2c(=O)[nH]1)[CH]")
    cytosine = Chem.MolFromSmarts("nc1cc[nH]c(=O)n1")
    thymine = Chem.MolFromSmarts("Cc1c[nH]c(=O)[nH]c1=O")
    uracil = Chem.MolFromSmarts("O=c1ccn(C)c(=O)[nH]1")

    # Check for nucleobase
    base_match = mol.HasSubstructMatch(adenine) or mol.HasSubstructMatch(guanine) or \
                 mol.HasSubstructMatch(cytosine) or mol.HasSubstructMatch(thymine) or \
                 mol.HasSubstructMatch(uracil)
    if not base_match:
        return False, "No nucleobase found"

    # Define sugar SMARTS patterns
    ribose = Chem.MolFromSmarts("[CR1][CR1]([OR1])[OR1][CR1]([OR1])[CR1]([OR1])[OR1]")
    deoxyribose = Chem.MolFromSmarts("[CR1][CR1]([OR1])[OR1][CR1]([OR1])[CR1]([HR1])[OR1]")
    cyclic_phosphate = Chem.MolFromSmarts("[CR1]1[OR1][CR1]([OR1])[CR1]([OR1])[CR1]([OP])[OR1]1")

    # Check for sugar moiety
    sugar_match = mol.HasSubstructMatch(ribose) or mol.HasSubstructMatch(deoxyribose) or \
                  mol.HasSubstructMatch(cyclic_phosphate)
    if not sugar_match:
        return False, "No sugar moiety found"

    # Check for phosphate group
    phosphate_pattern = Chem.MolFromSmarts("[PX4]([OX2])([OX1])([OX2])[OX2]")
    phosphate_match = mol.HasSubstructMatch(phosphate_pattern)
    if not phosphate_match:
        return False, "No phosphate group found"

    # Check for glycosidic bond between base and sugar
    glycosidic_bond_pattern = Chem.MolFromSmarts("[CR1][OR1][CR1]")
    glycosidic_bond_match = mol.HasSubstructMatch(glycosidic_bond_pattern)
    if not glycosidic_bond_match:
        return False, "No glycosidic bond found"

    return True, "Contains a nucleobase, sugar moiety, phosphate group, and glycosidic bond"