"""
Classifies: CHEBI:61384 sulfolipid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sulfolipid(smiles: str):
    """
    Determines if a molecule is a sulfolipid based on its SMILES string.
    A sulfolipid is a compound containing a sulfonic acid residue joined by a sulfur-oxygen bond (sulfate ester) to a lipid.

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

    # Look for sulfate ester group connected via oxygen-sulfur-oxygen-carbon bonds
    # SMARTS pattern for O-S(=O)(=O)-O-C
    sulfate_ester_pattern = Chem.MolFromSmarts("O[S](=O)(=O)O[C]")
    if not mol.HasSubstructMatch(sulfate_ester_pattern):
        return False, "No sulfate ester group connected to carbon found"

    # Check for long hydrocarbon chains to assess lipid characteristic
    # Look for aliphatic chains of at least 12 carbons (could be adjusted as needed)
    alkyl_chain_pattern = Chem.MolFromSmarts("C" + "C" * 11)  # Chain of 12 carbons
    if not mol.HasSubstructMatch(alkyl_chain_pattern):
        return False, "No long hydrocarbon chain found"

    return True, "Contains sulfate ester group connected via oxygen-sulfur bond and lipid characteristics"