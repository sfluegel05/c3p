"""
Classifies: CHEBI:32952 amine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_amine(smiles: str):
    """
    Determines if a molecule is an amine based on its SMILES string.
    An amine is a compound derived from ammonia by replacing one, two or three
    hydrogen atoms by hydrocarbyl groups

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for at least one nitrogen atom
    has_nitrogen = any(atom.GetAtomicNum() == 7 for atom in mol.GetAtoms())
    if not has_nitrogen:
        return False, "No nitrogen atom found"

    # General check for amine (N bonded to H or C)
    # [N;H0,H1,H2] is a nitrogen with 0, 1, or 2 Hs. [!#0] is a non-hydrogen atom.
    amine_pattern = Chem.MolFromSmarts("[N;H0,H1,H2][!#0]")
    if not mol.HasSubstructMatch(amine_pattern):
         return False, "No amine group found: Nitrogen not bonded to non-Hydrogen"

    # Exclude amides (-C(=O)-N-) using a substructure search to eliminate false positives
    amide_pattern = Chem.MolFromSmarts("[CX3](=[OX1])N")
    if mol.HasSubstructMatch(amide_pattern):
       return False, "Contains amide group"

    # If none of the above exclusion applies, it's an amine.
    return True, "Contains a nitrogen atom bonded to at least one non-hydrogen atom and is not an amide."