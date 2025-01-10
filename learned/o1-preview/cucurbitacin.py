"""
Classifies: CHEBI:16219 cucurbitacin
"""
"""
Classifies: cucurbitacin
"""

from rdkit import Chem
from rdkit.Chem import Descriptors

def is_cucurbitacin(smiles: str):
    """
    Determines if a molecule is a cucurbitacin based on its SMILES string.
    Cucurbitacins are tetracyclic triterpenoids derived from cucurbitane.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cucurbitacin, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for tetracyclic ring system
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings < 4:
        return False, f"Contains {num_rings} rings, less than 4 required for tetracyclic system"

    # Get the smallest set of smallest rings (SSSR)
    ssr = Chem.GetSymmSSSR(mol)
    if len(ssr) < 4:
        return False, "Less than 4 rings detected in the molecule"

    # Check ring sizes
    ring_sizes = [len(ring) for ring in ssr]
    num_6_rings = ring_sizes.count(6)
    num_5_rings = ring_sizes.count(5)
    if num_6_rings < 3 or num_5_rings < 1:
        return False, f"Incorrect ring sizes: {num_6_rings} six-membered rings and {num_5_rings} five-membered rings"

    # Attempt to match the cucurbitane skeleton (simplified)
    # Cucurbitane skeleton SMARTS pattern
    cucurbitane_smarts = 'C1CCC2C(C1)CCC3C(C2)CC=C4C3CCCC4'  # Simplified cucurbitane core
    cucurbitane_pattern = Chem.MolFromSmarts(cucurbitane_smarts)
    if not mol.HasSubstructMatch(cucurbitane_pattern):
        return False, "Does not contain the cucurbitane tetracyclic skeleton"

    # Check for characteristic functional groups
    # Hydroxyl groups: [OX2H]
    num_hydroxyl = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[OX2H]')))
    # Ketone groups: [CX3](=O)[#6]
    num_ketone = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[CX3](=O)[#6]')))
    if num_hydroxyl < 1 or num_ketone < 1:
        return False, "Does not contain the characteristic hydroxyl and ketone groups"

    # Check molecular weight
    mol_wt = Descriptors.ExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for a typical cucurbitacin"

    # Passed all checks
    return True, "Contains cucurbitane skeleton with characteristic functional groups"

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'cucurbitacin',
        'definition': 'Any one of a class of tetracyclic triterpenoids, formally derived from the triterpene hydrocarbon cucurbitane, developed by some plants (especially those of the family Cucurbitaceaeas) as a defence mechanism against herbivores.',
        'parents': []
    }
}