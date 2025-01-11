"""
Classifies: CHEBI:72588 semisynthetic derivative
"""
"""
Classifies: CHEBI:90773 semisynthetic derivative
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

def is_semisynthetic_derivative(smiles: str):
    """
    Determines if a molecule is a semisynthetic derivative based on its SMILES string.
    A semisynthetic derivative is an organic molecule derived from a natural product by partial chemical synthesis.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a semisynthetic derivative, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate Murcko scaffold of the molecule
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    if scaffold.GetNumAtoms() == 0:
        return False, "Unable to generate scaffold"

    # Define functional groups considered as synthetic modifications
    synthetic_groups = [
        Chem.MolFromSmarts('[#6][F,Cl,Br,I]'),  # Halogenated carbons
        Chem.MolFromSmarts('c[F,Cl,Br,I]'),     # Halogenated aromatics
        Chem.MolFromSmarts('N(=O)[O-]'),        # Nitro group
        Chem.MolFromSmarts('C#C'),              # Alkyne
        Chem.MolFromSmarts('N=[N+]=[N-]'),      # Azide
        Chem.MolFromSmarts('C(=O)O[C,N]'),      # Esters and carbamates
        Chem.MolFromSmarts('C(=O)N[C,N]'),      # Amides
        Chem.MolFromSmarts('CNS(=O)(=O)'),      # Sulfonamides
        Chem.MolFromSmarts('P(=O)(O)O'),        # Phosphate esters
        Chem.MolFromSmarts('C(=O)O[C@H]'),      # Esterification at chiral centers
        Chem.MolFromSmarts('[C;!R]=[C;!R]'),    # Unconjugated double bonds
    ]

    # Check for synthetic modifications
    has_synthetic_modification = False
    modifications = []

    for group in synthetic_groups:
        if mol.HasSubstructMatch(group):
            has_synthetic_modification = True
            smarts = Chem.MolToSmarts(group)
            modifications.append(smarts)

    # Compare scaffold to molecule
    if mol.GetNumAtoms() == scaffold.GetNumAtoms():
        return False, "Molecule appears to be a natural product without modifications"

    # Check for significant difference between scaffold and molecule
    diff = mol.GetNumAtoms() - scaffold.GetNumAtoms()
    if diff < 5 and not has_synthetic_modification:
        return False, "No significant synthetic modifications detected"

    # If synthetic modifications are present and scaffold is smaller than molecule, classify as semisynthetic
    if has_synthetic_modification:
        reason = f"Contains natural product scaffold with synthetic modifications: {', '.join(modifications)}"
        return True, reason
    else:
        return False, "No synthetic modifications detected"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:90773',
        'name': 'semisynthetic derivative',
        'definition': 'Any organic molecular entity derived from a natural product by partial chemical synthesis.',
        'parents': []
    }
}