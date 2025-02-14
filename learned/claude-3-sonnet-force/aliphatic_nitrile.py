"""
Classifies: CHEBI:80291 aliphatic nitrile
"""
"""
Classifies: CHEBI:33569 aliphatic nitrile
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_aliphatic_nitrile(smiles: str):
    """
    Determines if a molecule is an aliphatic nitrile based on its SMILES string.
    An aliphatic nitrile is a compound containing a nitrile (-C≡N) group and
    an aliphatic (non-aromatic) carbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an aliphatic nitrile, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a nitrile group
    nitrile_pattern = Chem.MolFromSmarts("C#N")
    has_nitrile = bool(mol.HasSubstructMatch(nitrile_pattern))

    if not has_nitrile:
        return False, "No nitrile group found"

    # Check for the presence of an aliphatic carbon chain
    aliphatic_chain_pattern = Chem.MolFromSmarts("[CH2,CH3]~[CH2,CH3]")
    has_aliphatic_chain = bool(mol.HasSubstructMatch(aliphatic_chain_pattern))

    if not has_aliphatic_chain:
        return False, "No aliphatic carbon chain found"

    # Check for common functional groups or substituents
    functional_groups = ["[OX1]", "[OX2H]", "[NX3]", "[NX3H]", "[SX2]", "[PX3]", "[PX4]"]
    has_functional_group = any(mol.HasSubstructMatch(Chem.MolFromSmarts(group)) for group in functional_groups)

    # Classify based on the presence of a nitrile group, aliphatic chain, and common functional groups
    if has_nitrile and has_aliphatic_chain and has_functional_group:
        return True, "Molecule contains a nitrile group and an aliphatic carbon chain with common functional groups"
    elif has_nitrile and has_aliphatic_chain:
        return True, "Molecule contains a nitrile group and an aliphatic carbon chain"
    else:
        return False, "Molecule does not meet the criteria for an aliphatic nitrile"