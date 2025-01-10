"""
Classifies: CHEBI:192499 anthoxanthin
"""
"""
Classifies: anthoxanthin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_anthoxanthin(smiles: str):
    """
    Determines if a molecule is an anthoxanthin based on its SMILES string.
    Anthoxanthins are a type of flavonoid pigments in plants, including flavones and flavonols.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an anthoxanthin, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for flavone and flavonol cores
    flavone_pattern = Chem.MolFromSmarts('[O]=C1C=CC(=CC1)c1ccccc1')  # 2-phenylchromen-4-one
    flavonol_pattern = Chem.MolFromSmarts('[O]=C1C=C(O)C=CC1c1ccccc1')  # 3-hydroxyflavone

    # Check for flavone core
    if mol.HasSubstructMatch(flavone_pattern):
        return True, "Contains flavone core structure"

    # Check for flavonol core
    if mol.HasSubstructMatch(flavonol_pattern):
        return True, "Contains flavonol core structure"

    # Check for glycosides of flavone or flavonol
    # Define SMARTS pattern for O-glycosylation at position 3 or 7
    glycosylation_pattern = Chem.MolFromSmarts('[O]-[*]')  # Simplified pattern for glycosidic bond
    flavone_oglycoside = Chem.CombineMols(flavone_pattern, glycosylation_pattern)
    flavonol_oglycoside = Chem.CombineMols(flavonol_pattern, glycosylation_pattern)

    if mol.HasSubstructMatch(flavone_oglycoside):
        return True, "Contains flavone O-glycoside structure"

    if mol.HasSubstructMatch(flavonol_oglycoside):
        return True, "Contains flavonol O-glycoside structure"

    # Check for C-glycosides at position 6 or 8
    # Define SMARTS pattern for C-glycosylation
    c_glycosylation_pattern = Chem.MolFromSmarts('[C]-[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O')

    # Positions 6 and 8 on flavone core
    flavone_cglycoside6 = Chem.MolFromSmarts('[O]=C1C=CC(=CC1-[C])c1ccccc1')
    flavone_cglycoside8 = Chem.MolFromSmarts('[O]=C1C=CC(=CC1C=C-[C])c1ccccc1')

    if mol.HasSubstructMatch(flavone_cglycoside6):
        if mol.HasSubstructMatch(c_glycosylation_pattern):
            return True, "Contains flavone C-glycoside at position 6"

    if mol.HasSubstructMatch(flavone_cglycoside8):
        if mol.HasSubstructMatch(c_glycosylation_pattern):
            return True, "Contains flavone C-glycoside at position 8"

    # Check for methylated derivatives
    # Define SMARTS pattern for methoxy groups attached to aromatic rings
    methoxy_pattern = Chem.MolFromSmarts('c-OC')

    if mol.HasSubstructMatch(flavone_pattern) or mol.HasSubstructMatch(flavonol_pattern):
        if mol.HasSubstructMatch(methoxy_pattern):
            return True, "Contains methylated flavone or flavonol core"

    # No anthoxanthin core structure found
    return False, "Does not contain anthoxanthin core structure"