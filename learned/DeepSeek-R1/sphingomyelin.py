"""
Classifies: CHEBI:64583 sphingomyelin
"""
"""
Classifies: CHEBI:17892 sphingomyelin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_sphingomyelin(smiles: str):
    """
    Determines if a molecule is a sphingomyelin based on its SMILES string.
    Sphingomyelins have a sphingoid base with an amide-linked fatty acid and a phosphorylcholine group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sphingomyelin, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for phosphorylcholine group
    pcho_pattern = Chem.MolFromSmarts("[O-]P(=O)(OCC[N+](C)(C)C)O[#6]")
    pcho_matches = mol.GetSubstructMatches(pcho_pattern)
    if len(pcho_matches) != 1:
        return False, f"Found {len(pcho_matches)} phosphorylcholine groups, need exactly 1"

    # Get sphingoid base connection point (carbon attached to phosphate oxygen)
    sphingoid_carbon_idx = pcho_matches[0][-1