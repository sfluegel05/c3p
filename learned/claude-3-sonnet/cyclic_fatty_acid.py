"""
Classifies: CHEBI:59238 cyclic fatty acid
"""
"""
Classifies: CHEBI:35831 cyclic fatty acid

A cyclic fatty acid is defined as any fatty acid containing anywhere in its structure a ring of atoms.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cyclic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a cyclic fatty acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cyclic fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (-C(=O)O)
    acid_pattern = Chem.MolFromSmarts("[C](=O)(O)")
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxylic acid group found"

    # Check for ring structure within the carbon chain
    chain_ring_pattern = Chem.MolFromSmarts("[R2]@[R2]@[R2]@[R2]")
    if not mol.HasSubstructMatch(chain_ring_pattern):
        return False, "No ring structure found within the carbon chain"

    # Check for long aliphatic carbon chain (fatty acid)
    aliphatic_chain_pattern = Chem.MolFromSmarts("[C;H3][C;H2]~[C;H2]~[C;H2]~[C;H2]~[C;H2]~[C;H2]~[C;H2]")
    if not mol.HasSubstructMatch(aliphatic_chain_pattern):
        return False, "No long aliphatic carbon chain found"

    # Check for unsaturation (optional)
    unsaturation_pattern = Chem.MolFromSmarts("[C]=[C]")
    unsaturation_count = len(mol.GetSubstructMatches(unsaturation_pattern))
    if unsaturation_count == 0:
        reason = "Saturated cyclic fatty acid"
    else:
        reason = f"Unsaturated cyclic fatty acid with {unsaturation_count} double bonds"

    return True, reason