"""
Classifies: CHEBI:61384 sulfolipid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sulfolipid(smiles: str):
    """
    Determines if a molecule is a sulfolipid based on its SMILES string.
    A sulfolipid is a compound containing a sulfonic acid residue joined by a carbon-sulfur bond to a lipid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sulfolipid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for sulfonic acid group connected via carbon-sulfur bond
    # Sulfonic acid group: S(=O)(=O)-OH attached to carbon via S-C bond
    # SMARTS pattern matching both protonated (OH) and deprotonated (O-) forms
    sulfonic_acid_pattern = Chem.MolFromSmarts("[#6][S](=O)(=O)[O;H1,-]")

    if not mol.HasSubstructMatch(sulfonic_acid_pattern):
        return False, "No sulfonic acid group connected via carbon-sulfur bond found"

    # Check for long hydrocarbon chain to assess lipid characteristic
    # Look for aliphatic chain of at least 12 continuous carbons
    # This pattern matches a chain of 12 or more non-ring sp3 carbons
    alkyl_chain_pattern = Chem.MolFromSmarts("[C;D3,D4;H2,H3][C;D2,D3;H2,H3]{10,}[C;D3,D4;H3]")  # Chain of at least 12 carbons

    if not mol.HasSubstructMatch(alkyl_chain_pattern):
        return False, "No long hydrocarbon chain of at least 12 carbons found"

    return True, "Contains sulfonic acid group connected via carbon-sulfur bond and lipid characteristics"