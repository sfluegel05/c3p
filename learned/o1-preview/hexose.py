"""
Classifies: CHEBI:18133 hexose
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_hexose(smiles: str):
    """
    Determines if a molecule is a hexose based on its SMILES string.
    A hexose is any six-carbon monosaccharide which in its linear form contains either
    an aldehyde group at position 1 (aldohexose) or a ketone group at position 2 (ketohexose).
    This function accounts for both linear and cyclic forms.

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

    # Count number of carbons
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons != 6:
        return False, f"Number of carbon atoms ({num_carbons}) is not 6"

    # Check for monosaccharide structure (no glycosidic bonds)
    # Count the number of rings to exclude disaccharides and polysaccharides
    if mol.GetRingInfo().NumRings() > 1:
        return False, "Molecule has more than one ring, may not be a monosaccharide"

    # Define SMARTS patterns for aldohexoses and ketohexoses in linear and cyclic forms
    # Aldohexose linear form
    aldohexose_linear = Chem.MolFromSmarts("[#6H]=O-[#6]-[#6]-[#6]-[#6]-[#6]")
    # Aldohexose cyclic form (hemiacetal pyranose and furanose)
    aldohexose_cyclic = Chem.MolFromSmarts("C1OC([#6H])(CO)C(O)C1O")
    # Ketohexose linear form
    ketohexose_linear = Chem.MolFromSmarts("[#6]-C(=O)-[#6]-[#6]-[#6]-[#6]-O")
    # Ketohexose cyclic form (hemiketal pyranose and furanose)
    ketohexose_cyclic = Chem.MolFromSmarts("C1OC(CO)C(O)C(O)C1O")

    # Match patterns
    if mol.HasSubstructMatch(aldohexose_linear):
        return True, "Molecule is an aldohexose in linear form"
    elif mol.HasSubstructMatch(ketohexose_linear):
        return True, "Molecule is a ketohexose in linear form"
    elif mol.HasSubstructMatch(aldohexose_cyclic):
        return True, "Molecule is an aldohexose in cyclic form"
    elif mol.HasSubstructMatch(ketohexose_cyclic):
        return True, "Molecule is a ketohexose in cyclic form"

    # As an alternative, check for general hexose pattern (both cyclic and linear)
    # Six carbons with hydroxyl groups and one oxygen in ring (for cyclic forms)
    hexose_pattern = Chem.MolFromSmarts("""
    [
        $([C;H1,H2][O;H1,H0]),
        $([C;H1,H2][O;H1,H0]),
        $([C;H1,H2][O;H1,H0]),
        $([C;H1,H2][O;H1,H0]),
        $([C;H1,H2][O;H1,H0]),
        $([C;H1,H2]=[O,N])
    ]
    """)

    if mol.HasSubstructMatch(hexose_pattern):
        return True, "Molecule matches general hexose pattern"

    return False, "Molecule does not match hexose patterns"