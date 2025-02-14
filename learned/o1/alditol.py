"""
Classifies: CHEBI:17522 alditol
"""
"""
Classifies: CHEBI:37277 alditol
"""
from rdkit import Chem

def is_alditol(smiles: str):
    """
    Determines if a molecule contains an alditol substructure based on its SMILES string.
    An alditol is an acyclic polyol having the general formula HOCH2[CH(OH)]nCH2OH,
    formally derivable from an aldose by reduction of the carbonyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains an alditol substructure, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the alditol pattern using SMARTS
    # Terminal CH2OH group: [C;H2;!R][O;H1]
    # Internal carbons: [C;H1,H2;!R]([O;H1]?)([O;H1]?)
    # Build patterns for chain lengths from 2 to 10 (adjustable range)
    min_chain_length = 2  # Minimum number of carbons in the chain
    max_chain_length = 10  # Maximum number of carbons in the chain

    alditol_found = False
    for chain_length in range(min_chain_length, max_chain_length + 1):
        # Start with terminal CH2OH group
        pattern = '[C;H2;!R]-[O;H1]'
        # Add internal carbons: allow for optional hydroxyl groups to account for deoxy sugars
        for _ in range(chain_length - 2):
            pattern += '-[C;H1,H2;!R]([O;H1]?)-[H]'
        # End with terminal CH2OH group
        pattern += '-[C;H2;!R]-[O;H1]'

        # Create the RDKit molecule from SMARTS
        alditol_pattern = Chem.MolFromSmarts(pattern)
        if alditol_pattern is None:
            continue  # Skip invalid patterns

        # Search for the pattern in the molecule
        matches = mol.GetSubstructMatches(alditol_pattern)
        if matches:
            alditol_found = True
            return True, f"Molecule contains an alditol substructure of chain length {chain_length}"

    if not alditol_found:
        return False, "No alditol substructure found"