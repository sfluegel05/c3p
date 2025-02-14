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
    # The pattern matches an acyclic chain of carbons where each carbon is sp3 hybridized,
    # connected to exactly one hydroxyl group, and terminal carbons are CH2OH groups.
    # We consider chains of varying lengths (n ≥ 1).

    # Generate patterns for chain lengths from 2 to 10 (adjustable range)
    min_chain_length = 2  # Minimum number of carbons in the chain (excluding terminal CH2OH)
    max_chain_length = 10  # Maximum number of carbons in the chain

    alditol_found = False
    for chain_length in range(min_chain_length, max_chain_length + 1):
        # Build SMARTS pattern
        # Start with terminal CH2OH group
        pattern = '[C;H2;X4;!R]-[O;X2;H1]'
        # Add internal CH(OH) groups
        pattern += ''.join(['-[C;H1;X4;!R]-[O;X2;H1]' for _ in range(chain_length - 2)])
        # Add ending CH2OH group
        pattern += '-[C;H2;X4;!R]-[O;X2;H1]'

        # Create the RDKit molecule from SMARTS
        alditol_pattern = Chem.MolFromSmarts(pattern)
        if alditol_pattern is None:
            continue  # Skip invalid patterns

        # Search for the pattern in the molecule
        if mol.HasSubstructMatch(alditol_pattern):
            alditol_found = True
            break  # Alditol substructure found

    if alditol_found:
        return True, f"Molecule contains an alditol substructure of chain length {chain_length}"
    else:
        return False, "No alditol substructure found"