"""
Classifies: CHEBI:61778 triterpenoid saponin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_triterpenoid_saponin(smiles: str):
    """
    Determines if a molecule is a triterpenoid saponin based on its SMILES string.
    A triterpenoid saponin is a terpene glycoside in which the terpene moiety is a triterpenoid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a triterpenoid saponin, False otherwise
        str: Reason for classification
    """

    # Parse SMILES into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # General pattern for hexopyranose sugars (glycosidic part)
    sugar_smarts = Chem.MolFromSmarts("O[C@@H]1[C@H](O)[C@@H](O)[C@@H](O)[C@H](O)[C@H]1O")
    if not mol.HasSubstructMatch(sugar_smarts):
        return False, "No hexopyranose sugar moiety found"

    # Generic pentacyclic triterpenoid backbone patterns
    triterpenoid_smarts_patterns = [
        Chem.MolFromSmarts("C1CC2CCCC3C4CCC5=CC(=O)CCC5C4CCC3C2C1"),  # Example pattern for typical triterpenoid
        Chem.MolFromSmarts("C1CCC2C(C1)C3CCC4C(C3C2)CCC5C4CCC6(C5)CC7=C(C=C6)C(=C7)C"),
    ]

    triterpenoid_found = False
    for pattern in triterpenoid_smarts_patterns:
        if mol.HasSubstructMatch(pattern):
            triterpenoid_found = True
            break

    if not triterpenoid_found:
        return False, "No triterpenoid backbone found"

    # Check for linkage between glycoside and triterpenoid (assuming ether linkage for simplicity)
    linkage_found = mol.HasSubstructMatch(Chem.MolFromSmarts("C-O-C"))
    if not linkage_found:
        return False, "No glycosidic linkage found"

    # Assess molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 600:
        return False, "Molecular weight too low for triterpenoid saponin"

    return True, "Structure matches a triterpenoid saponin with a glycosidic linkage"