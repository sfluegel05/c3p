"""
Classifies: CHEBI:53156 polychlorobiphenyl
"""
from rdkit import Chem

def is_polychlorobiphenyl(smiles: str):
    """
    Determines if a molecule is a polychlorinated biphenyl (PCB) based on its SMILES string.
    A PCB is a biphenyl compound containing between 2 and 10 chlorine atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a PCB, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for biphenyl structure (two connected benzene rings)
    # Use a more restrictive SMARTS to define a biphenyl core where all atoms in the biphenyl
    # are either C or H. This will avoid matching substructures that are not biphenyls.
    biphenyl_pattern = Chem.MolFromSmarts("[cH]1[cH][cH][cH][cH][cH]1-[cH]2[cH][cH][cH][cH][cH]2")
    if not mol.HasSubstructMatch(biphenyl_pattern):
        return False, "Not a biphenyl structure"

    # Make sure that there are no rings connected to the biphenyl
    # using a single bond
    ring_pattern = Chem.MolFromSmarts("[cH]1[cH][cH][cH][cH][cH]1-[!#6,#1!]")
    if mol.HasSubstructMatch(ring_pattern):
        return False, "Other rings directly connected to biphenyl core"

    # Count chlorine atoms
    chlorine_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 17)

    # Check if chlorine count is within the range 2-10
    if chlorine_count < 2:
         return False, f"Too few chlorine atoms: {chlorine_count}. Must be between 2 and 10."
    if chlorine_count > 10:
        return False, f"Too many chlorine atoms: {chlorine_count}. Must be between 2 and 10."

    return True, "Biphenyl structure with 2-10 chlorine atoms"