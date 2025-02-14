"""
Classifies: CHEBI:36498 galactosylceramide
"""
"""
Classifies: CHEBI:34219 galactosylceramide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_galactosylceramide(smiles: str):
    """
    Determines if a molecule is a galactosylceramide based on its SMILES string.
    A galactosylceramide is a ceramide with a galactose monosaccharide head group.
    It consists of a sphingoid base amide-linked to a fatty acid, with a galactose
    connected via a glycosidic bond to the primary hydroxyl of the sphingoid base.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a galactosylceramide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define sphingoid base pattern (sphingosine, sphinganine, phytosphingosine)
    sphingoid_smarts = """
    [#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#8]  # long chain with terminal hydroxyl
    """
    sphingoid_pattern = Chem.MolFromSmarts("""
    [#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#7]-[#6]-[#8]
    """)

    # Generalized sphingoid base pattern: long-chain amino alcohol (allowing variable chain lengths and unsaturation)
    sphingoid_pattern = Chem.MolFromSmarts("""
    [#6](-[#6]){12,}(-[#6])-[#6]-[#6]-[#7]-[#6]-[#8]
    """)

    # Identify sphingoid base
    if not mol.HasSubstructMatch(sphingoid_pattern):
        return False, "No sphingoid base found"

    # Identify amide bond to fatty acid
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide bond to fatty acid found"

    # Identify galactose head group (both alpha and beta forms)
    galactose_pattern = Chem.MolFromSmarts("""
    [C@@H]1([C@@H]([C@H]([C@@H](CO1)O)O)O)O  # beta-D-galactose
    """)
    galactose_pattern_alpha = Chem.MolFromSmarts("""
    [C@H]1([C@@H]([C@H]([C@H](CO1)O)O)O)O  # alpha-D-galactose
    """)
    has_galactose = mol.HasSubstructMatch(galactose_pattern) or mol.HasSubstructMatch(galactose_pattern_alpha)
    if not has_galactose:
        return False, "No galactose head group found"

    # Check for glycosidic bond between sphingoid base and galactose
    glycosidic_bond_pattern = Chem.MolFromSmarts("[C@H,@@H](CO[C@H,@@H]1O[C@H,@@H]([C@H,@@H]([C@H,@@H](C1)O)O)O)[NH]")
    if not mol.HasSubstructMatch(glycosidic_bond_pattern):
        return False, "No glycosidic bond between sphingoid base and galactose found"

    return True, "Molecule is a galactosylceramide with sphingoid base amide-linked to fatty acid and galactose head group"