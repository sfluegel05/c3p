"""
Classifies: CHEBI:61655 steroid saponin
"""
"""
Classifies: CHEBI:37668 steroid saponin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_steroid_saponin(smiles: str):
    """
    Determines if a molecule is a steroid saponin based on its SMILES string.
    A steroid saponin is a molecule with a steroid backbone and one or more sugar moieties.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a steroid saponin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more general steroid backbone pattern (tetracyclic structure)
    steroid_pattern = Chem.MolFromSmarts("[C@]12[C@@]3([C@@]([C@@H]4[C@@]([C@@H](CC4)C)(C)C)(CC3)C)[C@@H](C1)CC2")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"

    # Define a more flexible sugar moiety pattern (glycosidic bond)
    sugar_pattern = Chem.MolFromSmarts("[C@@H]1([C@@H]([C@H]([C@@H]([C@H](O1)O)O)O)O)")
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if len(sugar_matches) == 0:
        return False, "No sugar moieties found"

    # Check for glycosidic bonds (O-C bond between steroid and sugar)
    glycosidic_bond_pattern = Chem.MolFromSmarts("[C@]12[C@@]3([C@@]([C@@H]4[C@@]([C@@H](CC4)C)(C)C)(CC3)C)[C@@H](C1)CC2.O[C@@H]1[C@@H]([C@H]([C@@H]([C@H](O1)O)O)O)O")
    if not mol.HasSubstructMatch(glycosidic_bond_pattern):
        return False, "No glycosidic bond found between steroid and sugar"

    # Count the number of sugar moieties
    n_sugars = len(sugar_matches)
    if n_sugars < 1:
        return False, "At least one sugar moiety is required"

    # Check molecular weight - steroid saponins typically have high molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, "Molecular weight too low for steroid saponin"

    return True, "Contains steroid backbone with one or more sugar moieties attached via glycosidic bonds"