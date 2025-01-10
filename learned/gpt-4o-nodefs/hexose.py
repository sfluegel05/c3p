"""
Classifies: CHEBI:18133 hexose
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

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

    # Ensure molecule is organic
    if not all(atom.GetAtomicNum() in (1, 6, 7, 8, 16) for atom in mol.GetAtoms()):
        return False, "Non-organic elements present"

    # Check if the molecule is an ion or salt and reduce it to the main organic component
    fragments = Chem.GetMolFrags(mol, asMols=True)
    main_mol = max(fragments, key=lambda frag: frag.GetNumAtoms())
    
    # Check for 6-carbon backbone (especially for hexoses)
    num_carbons = sum(1 for atom in main_mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons != 6:
        return False, f"Number of carbon atoms is {num_carbons}, expected 6 for hexose"

    # Check for presence of hydroxyl groups (typically in hexoses)
    num_hydroxyls = sum(1 for atom in main_mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetDegree() == 1)
    if num_hydroxyls < 4:
        return False, f"Insufficient hydroxyl groups, found {num_hydroxyls}, expected at least 4"

    # Identify cyclic structures for pyranoses and furanoses
    pyranose_pattern = Chem.MolFromSmarts("C1(O)COC(OC1)[O]")  # broader pyranose
    furanose_pattern = Chem.MolFromSmarts("C1(O)COC(O)C1")    # broader furanose
    if main_mol.HasSubstructMatch(pyranose_pattern) or main_mol.HasSubstructMatch(furanose_pattern):
        return True, "Hexose in cyclic pyranose or furanose form"

    # For open-chain forms, check presence of aldehyde or ketone
    aldehyde_or_ketone_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")  # indicative of aldose/ketose
    if main_mol.HasSubstructMatch(aldehyde_or_ketone_pattern):
        return True, "Hexose in open-chain form with carbonyl"

    # Check for N-acetylation or similar modifications
    n_acetyl_pattern = Chem.MolFromSmarts("NC(=O)C")  # common N-acetylation
    if main_mol.HasSubstructMatch(n_acetyl_pattern):
        return True, "Hexose derivative with N-acetyl modification"

    # Allow broader acetylated structures
    acetyl_patterns = [
        Chem.MolFromSmarts("OC(=O)C"),  # O-acetyl
        Chem.MolFromSmarts("C(=O)OC"),  # ester linkage
    ]
    for pattern in acetyl_patterns:
        if main_mol.HasSubstructMatch(pattern):
            return True, "Hexose derivative with acetylation"

    return False, "Structure does not fit typical or expanded hexose patterns"