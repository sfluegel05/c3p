"""
Classifies: CHEBI:18241 2'-deoxyribonucleoside 5'-monophosphate
"""
"""
Classifies: CHEBI:28885 2'-deoxyribonucleoside 5'-monophosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_2__deoxyribonucleoside_5__monophosphate(smiles: str):
    """
    Determines if a molecule is a 2'-deoxyribonucleoside 5'-monophosphate based on its SMILES string.
    A 2'-deoxyribonucleoside 5'-monophosphate has a 2'-deoxyribose sugar, a nucleobase, and a phosphate group at the 5' position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2'-deoxyribonucleoside 5'-monophosphate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a phosphate group (P=O or P-O, including deprotonated forms)
    phosphate_pattern = Chem.MolFromSmarts("[PX4](=O)([OX2,-])([OX2,-])")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) == 0:
        return False, "No phosphate group found"

    # Check for the presence of a 2'-deoxyribose sugar (C1C(C(COP(=O)(O)O)O)O)
    deoxyribose_pattern = Chem.MolFromSmarts("C1C(C(COP(=O)([OX2,-])[OX2,-])O)O1")
    deoxyribose_matches = mol.GetSubstructMatches(deoxyribose_pattern)
    if len(deoxyribose_matches) == 0:
        return False, "No 2'-deoxyribose sugar found"

    # Check for the presence of a nucleobase (adenine, guanine, cytosine, thymine, or uracil, including modifications)
    nucleobase_patterns = [
        Chem.MolFromSmarts("Nc1ncnc2n(cnc12)"),  # Adenine
        Chem.MolFromSmarts("Nc1nc2c(n1)[nH]c(=O)n2"),  # Guanine
        Chem.MolFromSmarts("Nc1ccnc(=O)n1"),  # Cytosine
        Chem.MolFromSmarts("Cc1cnc(=O)[nH]c1=O"),  # Thymine
        Chem.MolFromSmarts("O=C1NC(=O)C=CN1"),  # Uracil
        Chem.MolFromSmarts("Nc1nc2c(n1)[nH]c(=O)n2"),  # Modified Guanine
        Chem.MolFromSmarts("Nc1ccnc(=O)n1"),  # Modified Cytosine
        Chem.MolFromSmarts("Cc1cnc(=O)[nH]c1=O"),  # Modified Thymine
        Chem.MolFromSmarts("O=C1NC(=O)C=CN1")  # Modified Uracil
    ]
    nucleobase_found = any(mol.HasSubstructMatch(pattern) for pattern in nucleobase_patterns)
    if not nucleobase_found:
        return False, "No nucleobase found"

    # Check that the phosphate group is attached to the 5' position of the sugar
    phosphate_atom = phosphate_matches[0][0]
    sugar_atom = deoxyribose_matches[0][0]
    if not mol.GetBondBetweenAtoms(phosphate_atom, sugar_atom):
        return False, "Phosphate group not attached to the 5' position of the sugar"

    return True, "Contains a 2'-deoxyribose sugar, a nucleobase, and a phosphate group at the 5' position"