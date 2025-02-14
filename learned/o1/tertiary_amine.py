"""
Classifies: CHEBI:32876 tertiary amine
"""
"""
Classifies: CHEBI:32876 tertiary amine
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_tertiary_amine(smiles: str):
    """
    Determines if a molecule is a tertiary amine based on its SMILES string.
    A tertiary amine is a compound where a nitrogen atom is sp3-hybridized,
    bonded to three hydrocarbyl groups (carbon atoms), not part of an amide, imine,
    nitrile, nitro group, or aromatic system.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tertiary amine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for tertiary amine
    # Nitrogen atom that is:
    # - sp3 hybridized (implicit in [NX3])
    # - Not aromatic
    # - Has no hydrogen atoms attached (H0)
    # - Bonded to three carbon atoms
    # - Excludes nitrogen atoms in amides, imines, nitriles, nitro groups
    # - Neighboring carbons are not double-bonded to heteroatoms (e.g., carbonyl groups)
    tertiary_amine_smarts = "[NX3;H0;!$(N-C=O);!$(N-C=N);!$(N#C);!$(N[N+]=O);!a]([C;!$(C=[O,N,S])])([C;!$(C=[O,N,S])])[C;!$(C=[O,N,S])]"

    tertiary_amine_pattern = Chem.MolFromSmarts(tertiary_amine_smarts)

    # Check if the molecule matches the tertiary amine pattern
    if mol.HasSubstructMatch(tertiary_amine_pattern):
        return True, "Contains a tertiary amine group"
    else:
        return False, "No tertiary amine group found"