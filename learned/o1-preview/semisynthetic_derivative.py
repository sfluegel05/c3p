"""
Classifies: CHEBI:72588 semisynthetic derivative
"""
"""
Classifies: CHEBI:90773 semisynthetic derivative
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors

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

    # Define SMARTS patterns for common natural product scaffolds
    patterns = {
        'Steroid': Chem.MolFromSmarts('C1CCC2C(C1)CCC3C2CCC4(C3CCC4)'),
        'Beta-lactam': Chem.MolFromSmarts('C1[C@H]([C@@H](N1)C(=O)O)C(=O)O'),
        'Macrolide': Chem.MolFromSmarts('C1CC[C@H]2[C@@H](C1)O[C@H](C(=O))[C@H](O)[C@H](C)[C@@H]2O'),
        'Glycoside': Chem.MolFromSmarts('O[C@H]1[C@@H](O)[C@H](O)[C@H](O)[C@H]1O'),
        'Alkaloid': Chem.MolFromSmarts('C1CNCCN1'),
        'Polyketide': Chem.MolFromSmarts('C(=O)C(CC(=O))CC(=O)'),
        'Terpene': Chem.MolFromSmarts('C=C(C)CCC=C(C)C'),
        'Peptide': Chem.MolFromSmarts('N[C@@H](C)C(=O)N'),
    }

    # Check for natural product scaffolds
    has_scaffold = False
    scaffold_names = []
    for name, pattern in patterns.items():
        if mol.HasSubstructMatch(pattern):
            has_scaffold = True
            scaffold_names.append(name)

    if not has_scaffold:
        return False, "No natural product scaffold found"

    # Check for synthetic modifications
    # Look for presence of halogens (Cl, Br, F, I)
    halogen_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() in [9, 17, 35, 53]]
    total_halogens = len(halogen_atoms)

    # Look for nitro groups
    nitro_group = Chem.MolFromSmarts('[N+](=O)[O-]')
    has_nitro = mol.HasSubstructMatch(nitro_group)

    # Look for azide groups
    azide_group = Chem.MolFromSmarts('N=[N+]=[N-]')
    has_azide = mol.HasSubstructMatch(azide_group)

    # Look for alkyne groups (rare in natural products)
    alkyne_group = Chem.MolFromSmarts('C#C')
    has_alkyne = mol.HasSubstructMatch(alkyne_group)

    # Consider presence of unnatural functional groups
    synthetic_modifications = total_halogens + has_nitro + has_azide + has_alkyne

    # Determine if synthetic modifications are present
    if synthetic_modifications > 0:
        scaffold_list = ', '.join(scaffold_names)
        reasons = []
        if total_halogens > 0:
            halogens = ', '.join([atom.GetSymbol() for atom in halogen_atoms])
            reasons.append(f"Contains halogens ({halogens})")
        if has_nitro:
            reasons.append("Contains nitro group")
        if has_azide:
            reasons.append("Contains azide group")
        if has_alkyne:
            reasons.append("Contains alkyne group")
        reason_str = "; ".join(reasons)
        return True, f"Contains natural product scaffold(s) ({scaffold_list}) with synthetic modifications: {reason_str}"
    else:
        return False, "Contains natural product scaffold but no synthetic modifications detected"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:90773',
        'name': 'semisynthetic derivative',
        'definition': 'Any organic molecular entity derived from a natural product by partial chemical synthesis.',
        'parents': []
    }
}