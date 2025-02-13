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

    # Aromaticity perception (ensure correct handling of aromatic rings)
    Chem.SanitizeMol(mol)

    # Define SMARTS patterns for flavone and flavonol cores
    # Flavone core: 2-phenylchromen-4-one
    flavone_pattern = Chem.MolFromSmarts('c1cc(-c2oc3ccccc3c(=O)c2)ccc1')

    # Flavonol core: 3-hydroxyflavone
    flavonol_pattern = Chem.MolFromSmarts('c1cc(-c2oc3ccccc3c(=O)c2O)ccc1')

    # Check for flavone core
    if mol.HasSubstructMatch(flavone_pattern):
        return True, "Contains flavone core structure"

    # Check for flavonol core
    if mol.HasSubstructMatch(flavonol_pattern):
        return True, "Contains flavonol core structure"

    # Define patterns for glycosylated flavones/flavonols (O-glycosides)
    # Simplified pattern: flavone or flavonol core with any sugar attached via oxygen
    o_glycoside_pattern = Chem.MolFromSmarts('[$([cO][CX4H]),$([cO][CX4H][CX4H]),$([cO][CX4H][CX4H][CX4H])]')
    flavone_oglycoside_pattern = Chem.CombineMols(flavone_pattern, o_glycoside_pattern)
    flavonol_oglycoside_pattern = Chem.CombineMols(flavonol_pattern, o_glycoside_pattern)

    if mol.HasSubstructMatch(flavone_oglycoside_pattern):
        return True, "Contains flavone O-glycoside structure"

    if mol.HasSubstructMatch(flavonol_oglycoside_pattern):
        return True, "Contains flavonol O-glycoside structure"

    # Check for C-glycosides (sugar attached directly to carbon)
    # C-glycoside pattern: flavone or flavonol core with sugar attached to C6 or C8
    c_glycoside_pattern = Chem.MolFromSmarts('c1cc(-c2oc3ccccc3c(=O)c2[c,C])ccc1')
    if mol.HasSubstructMatch(c_glycoside_pattern):
        return True, "Contains flavone or flavonol C-glycoside structure"

    # Check for methylated derivatives (methoxy groups)
    methoxy_pattern = Chem.MolFromSmarts('c-oc')
    if mol.HasSubstructMatch(flavone_pattern) or mol.HasSubstructMatch(flavonol_pattern):
        if mol.HasSubstructMatch(methoxy_pattern):
            return True, "Contains methylated flavone or flavonol core"

    # Check for sulfonated derivatives (sulfo groups)
    sulfo_pattern = Chem.MolFromSmarts('S(=O)(=O)[O-]')
    if mol.HasSubstructMatch(sulfo_pattern):
        if mol.HasSubstructMatch(flavone_pattern) or mol.HasSubstructMatch(flavonol_pattern):
            return True, "Contains sulfonated flavone or flavonol"

    # No anthoxanthin core structure found
    return False, "Does not contain anthoxanthin core structure"