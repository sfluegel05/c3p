"""
Classifies: CHEBI:39437 tocol
"""
"""
Classifies: CHEBI:33235 tocol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_tocol(smiles: str):
    """
    Determines if a molecule is a tocol based on its SMILES string.
    A tocol is a chromanol with a chroman-6-ol skeleton that is substituted at position 2 by a saturated or triply-unsaturated hydrocarbon chain consisting of three isoprenoid units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tocol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for chroman-6-ol skeleton pattern
    # The pattern should match the chroman-6-ol skeleton with a hydroxyl group at position 6
    chromanol_pattern = Chem.MolFromSmarts("[O]1[C@@]([C])([C])([C])[C]2=C1C(=C([C])C(=C2)O)")
    if not mol.HasSubstructMatch(chromanol_pattern):
        return False, "No chroman-6-ol skeleton found"

    # Look for the hydrocarbon chain at position 2
    # The chain should consist of three isoprenoid units (C5H8) which can be saturated or triply unsaturated
    # We will search for a pattern that matches the general structure of the chain
    hydrocarbon_chain_pattern = Chem.MolFromSmarts("[C]~[C]~[C]~[C]~[C]~[C]~[C]~[C]~[C]~[C]~[C]~[C]~[C]~[C]~[C]")
    hydrocarbon_chain_matches = mol.GetSubstructMatches(hydrocarbon_chain_pattern)
    if len(hydrocarbon_chain_matches) == 0:
        return False, "No hydrocarbon chain found at position 2"

    # Check if the chain is connected to the chroman-6-ol skeleton at position 2
    # We will use a more specific pattern to ensure the chain is attached at the correct position
    tocol_pattern = Chem.MolFromSmarts("[O]1[C@@]([C])([C])([C])[C]2=C1C(=C([C])C(=C2)O)[C]~[C]~[C]~[C]~[C]~[C]~[C]~[C]~[C]~[C]~[C]~[C]~[C]~[C]~[C]")
    if not mol.HasSubstructMatch(tocol_pattern):
        return False, "Hydrocarbon chain not connected at position 2 of the chroman-6-ol skeleton"

    # Check the number of rotatable bonds to ensure the chain is long enough
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, "Hydrocarbon chain too short to be three isoprenoid units"

    # Check molecular weight - tocols typically have a molecular weight >400 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, "Molecular weight too low for tocol"

    return True, "Contains chroman-6-ol skeleton with a hydrocarbon chain of three isoprenoid units at position 2"