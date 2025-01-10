"""
Classifies: CHEBI:26377 pterocarpans
"""
"""
Classifies: CHEBI:26395 pterocarpans
"""
from rdkit import Chem

def is_pterocarpans(smiles: str):
    """
    Determines if a molecule is a pterocarpan based on its SMILES string.
    A pterocarpan is characterized by a 6a,11a-dihydro-6H-[1]benzofuro[3,2-c]chromene skeleton,
    which is a fused tricyclic system involving benzofuran and chromene rings.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a pterocarpan, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the pterocarpan core structure using a generalized SMARTS pattern
    # This pattern represents the fused benzofurochromene core
    pterocarpan_core_smarts = """
    [#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1
    [C@H]2COc3ccccc3O[C@@H]2
    """  # SMARTS for pterocarpan core without specifying substituents

    # Remove whitespace and newlines from SMARTS
    pterocarpan_core_smarts = ''.join(pterocarpan_core_smarts.split())

    # Convert SMARTS to molecule
    pterocarpan_core_mol = Chem.MolFromSmarts(pterocarpan_core_smarts)
    if pterocarpan_core_mol is None:
        return False, "Could not define pterocarpan core SMARTS"

    # Check for substructure match without considering chirality
    if mol.HasSubstructMatch(pterocarpan_core_mol, useChirality=False):
        return True, "Molecule contains the pterocarpan core structure"

    return False, "Molecule does not contain the pterocarpan core structure"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:26395',
        'name': 'pterocarpans',
        'definition': 'Members of the class of benzofurochromene with a 6a,11a-dihydro-6H-[1]benzofuro[3,2-c]chromene skeleton and its substituted derivatives. They generally bear structural resemblance to isoflavanoids that possess antibiotic activity and are produced by plant tissues in response to infection. They are the 3,4-dihydroderivatives of coumestans.',
        'parents': []
    }
}