"""
Classifies: CHEBI:18133 hexose
"""
from rdkit import Chem

def is_hexose(smiles: str):
    """
    Determines if a molecule is a hexose or a derivative thereof based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hexose, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule is an ion or salt and reduce it to the main organic component
    fragments = Chem.GetMolFrags(mol, asMols=True)
    main_mol = max(fragments, key=lambda frag: frag.GetNumAtoms())

    # Check for 6-carbon backbone (specifically focusing on hexoses)
    num_carbons = sum(1 for atom in main_mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons != 6:
        return False, f"Number of carbon atoms is {num_carbons}, expected 6 for hexose"

    # Check for presence of multiple hydroxyl groups (focus on typical hexose functionality)
    num_hydroxyls = sum(1 for atom in main_mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetDegree() == 1)
    if num_hydroxyls < 4:
        return False, f"Insufficient hydroxyl groups, found {num_hydroxyls}, expected at least 4"

    # Identify cyclic structures: pyranose or furanose
    pyranose_pattern = Chem.MolFromSmarts("C1(CO)OC(O)C(O)C(O)O1")  # pyranose
    furanose_pattern = Chem.MolFromSmarts("C1(CO)OC(O)C(O)C1O")    # furanose
    if main_mol.HasSubstructMatch(pyranose_pattern) or main_mol.HasSubstructMatch(furanose_pattern):
        return True, "Hexose in cyclic pyranose or furanose form"

    # For open-chain forms, ensure presence of aldehyde or ketone groups with hydroxyls
    aldehyde_or_ketone_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")  # =O group indicative of aldose/ketose
    if main_mol.HasSubstructMatch(aldehyde_or_ketone_pattern):
        return True, "Hexose in open-chain form with carbonyl"

    # Check for possible acetylation - no more than one acetyl group
    acetyl_pattern = Chem.MolFromSmarts("CC(=O)O")  # Common acetylation form
    num_acetyls = len(main_mol.GetSubstructMatches(acetyl_pattern))
    if num_acetyls > 1:
        return False, "Too many acetyl groups, structure diverges from hexose"

    return False, "Structure does not fit typical or derivative hexose patterns"