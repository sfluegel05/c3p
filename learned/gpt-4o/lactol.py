"""
Classifies: CHEBI:38131 lactol
"""
from rdkit import Chem

def is_lactol(smiles: str):
    """
    Determines if a molecule is a lactol based on its SMILES string.
    A lactol is a cyclic hemiacetal formed by the intramolecular addition
    of a hydroxyl group to an aldehydic or ketonic carbonyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lactol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define an improved lactol pattern:
    # The hemiacetal cyclic structure: an oxygen (as part of the ring, O1) bonded consecutively to:
    # C in ring (for hemiacetal formation), next to one C having [OH]
    lactol_pattern = Chem.MolFromSmarts('[#6][C@H]([O])[O]1[#6][#6]1') 

    # Match the lactol pattern
    if mol.HasSubstructMatch(lactol_pattern):
        return True, "Contains cyclic hemiacetal (lactol structure)"

    return False, "Lacks structural features of a lactol: cyclic hemiacetal"

# Test on identified lactol examples and review outcomes
examples = [
    "OC[C@]1(O)OC[C@@H](O)[C@@H](O)[C@@H]1O",  # alpha-D-fructopyranose
    "OC[C@@]1(O)OC[C@@H](O)[C@@H](O)[C@@H]1O",  # beta-D-fructopyranose
    "O1C2=C(C(CC1(C3=CC(=C(C=C3)O)O)O)=O)C(=CC(=C2)O)O",  # 2-(3,4-dihydroxyphenyl)-2,5,7-trihydroxy-2,3-dihydro-4H-chromen-4-one
]

for smiles in examples:
    result, reason = is_lactol(smiles)
    print(f"SMILES: {smiles}\nIs Lactol: {result}\nReason: {reason}\n")