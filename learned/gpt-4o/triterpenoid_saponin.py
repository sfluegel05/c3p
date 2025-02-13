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

    # Expand the glycoside detection: generic hexopyranose pattern
    glycosidic_smarts = Chem.MolFromSmarts("O[C@H]1[C@H](O)[C@H](O)[C@H](O)[C@@H](O)C1")
    if not mol.HasSubstructMatch(glycosidic_smarts):
        return False, "No applicable sugar moiety (glycoside) found"

    # Example of SMARTS pattern for a more versatile triterpenoid detection
    triterpenoid_smarts_patterns = [
        Chem.MolFromSmarts("C1CCC2C1CCC3C2CCC4C3CCC5C4CCCC5(C)C"),
        Chem.MolFromSmarts("C1CCC2C(C1)CCC3C2CCC4C3CCC5C4CCCC5C")
    ]

    triterpenoid_found = False
    for pattern in triterpenoid_smarts_patterns:
        if mol.HasSubstructMatch(pattern):
            triterpenoid_found = True
            break

    if not triterpenoid_found:
        return False, "No triterpenoid moiety found"

    # Check for an ester or ether linkage
    linkage_found = any(
        mol.HasSubstructMatch(Chem.MolFromSmarts(smarts))
        for smarts in ["O-C=O", "O-[C]-O"]
    )
    if not linkage_found:
        return False, "No linkage between glycoside and triterpenoid"

    # Assess molecule complexity
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 800:
        return False, "Molecular weight too low for typical triterpenoid saponin"

    return True, "Contains features of a triterpenoid saponin with glycoside linkage"