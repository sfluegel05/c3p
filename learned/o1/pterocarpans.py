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
    A pterocarpan is characterized by a 6a,11a-dihydro-6H-[1]benzofuro[3,2-c]chromene skeleton.

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

    # Define the pterocarpan core structure as a molecule
    pterocarpan_core_smiles = 'C1=CC=C2C(=C1)C3C(O2)COC4=C3C=CC=C4'  # SMILES for pterocarpan core
    pterocarpan_core_mol = Chem.MolFromSmiles(pterocarpan_core_smiles)
    if pterocarpan_core_mol is None:
        return False, "Could not define pterocarpan core structure"

    # Convert core structure to SMARTS pattern
    pterocarpan_core_smarts = Chem.MolToSmarts(pterocarpan_core_mol)

    # Check for substructure match
    if mol.HasSubstructMatch(Chem.MolFromSmarts(pterocarpan_core_smarts)):
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