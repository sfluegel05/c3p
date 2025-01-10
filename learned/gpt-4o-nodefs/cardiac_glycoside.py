"""
Classifies: CHEBI:83970 cardiac glycoside
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cardiac_glycoside(smiles: str):
    """
    Determines if a molecule is a cardiac glycoside based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule is a cardiac glycoside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Steroid core recognition (cyclopenta[a]phenanthrene skeleton)
    steroid_patterns = [
        Chem.MolFromSmarts('[C@@]12[C@@H]3CC[C@H]4[C@H]([C@](C3)(CC[C@@H]4[C@@H]1CC2)O)[C@H](C)C'),
        Chem.MolFromSmarts('[C@@H]1[C@@H]2C[C@@H]3CC[C@@H]4[C@]3(C2CC[C@@H]4C1)C')
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in steroid_patterns):
        return False, "No steroid core found"

    # Glycoside detection such as O-glycosidic bond
    glycoside_patterns = [
        Chem.MolFromSmarts('[CX4H]([OX2][#6])-[#6]-[#8]'),
        Chem.MolFromSmarts('O[C@H]1[C@H]([C@@H]([C@H]([C@H](CO1)O)O)O)O') # representative of sugar motifs
    ]
    if not any(mol.HasSubstructMatch(glycoside_pattern) for glycoside_pattern in glycoside_patterns):
        return False, "No glycosidic linkage found"

    # Lactone ring detection
    lactone_patterns = [
        Chem.MolFromSmarts('O=C1CCOC1'),
        Chem.MolFromSmarts('O=C1CCC=C1')  # Including more variants
    ]
    if not any(mol.HasSubstructMatch(lactone_pattern) for lactone_pattern in lactone_patterns) and \
       not any(mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)) for pattern in [
           '[C@]4([C@H](O)C(=O)[C@@H](OC)C=O)'])
        return False, "No lactone or related ester ring found"
    
    # Hydroxyl groups verification
    hydroxyl_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'O' and atom.GetTotalDegree() == 1)
    if hydroxyl_count < 3:
        return False, f"Insufficient hydroxyl groups in the molecule: {hydroxyl_count} found"
    
    return True, "The structure contains a steroid core with glycosidic linkage and a lactone ring characteristic of cardiac glycosides"