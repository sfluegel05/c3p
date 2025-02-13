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

    # General hexopyranose structure (more flexible than before)
    glycosidic_smarts = Chem.MolFromSmarts("C1(CO)C(O)C(O)C(O)C(O)O1")
    if not mol.HasSubstructMatch(glycosidic_smarts):
        return False, "No applicable sugar moiety (glycoside) found"

    # Expanded SMARTS patterns to match triterpenoid backbones
    triterpenoid_smarts_patterns = [
        Chem.MolFromSmarts("C1(CC2)(C3)(C4)C5CC(C)(C)C6CC(C)(C)CC7(C)C1CC1=CC[C@H](O[C@]1([H])O[C@@H](CC(C█))O[H])C8(C)C"),
        Chem.MolFromSmarts("C1CCCCC1C2CCCC(C)(C)C3CCCCC4(CC)CCCC5(C)CCCC(C4)C5"),
        # More generic pentacyclic terpenic scaffold
        Chem.MolFromSmarts("C1CC2CCCC3C2CCC4C3CCC5C4CCCC5(C)C"),
    ]

    triterpenoid_found = False
    for pattern in triterpenoid_smarts_patterns:
        if mol.HasSubstructMatch(pattern):
            triterpenoid_found = True
            break

    if not triterpenoid_found:
        return False, "No triterpenoid moiety found"

    # Check for a linkage between the glycoside and triterpenoid, covering esters or ethers
    linkage_found = any(
        mol.HasSubstructMatch(Chem.MolFromSmarts(smarts))
        for smarts in ["C(O)C(=O)", "C-O-C", "C(=O)O"]
    )
    if not linkage_found:
        return False, "No linkage between glycoside and triterpenoid"

    # Assess molecular complexity
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 600:
        return False, "Molecular weight too low for typical triterpenoid saponin"

    return True, "Contains features of a triterpenoid saponin with glycoside linkage"