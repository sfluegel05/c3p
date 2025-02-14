"""
Classifies: CHEBI:16412 lipopolysaccharide
"""
"""
Classifies: lipopolysaccharide
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_lipopolysaccharide(smiles: str):
    """
    Determines if a molecule is a lipopolysaccharide based on its SMILES string.
    A lipopolysaccharide consists of a lipid A moiety (N-acylated glucosamine disaccharide with fatty acids),
    a core polysaccharide (including KDO and heptoses), and an O-antigen polysaccharide.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lipopolysaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define pattern for N-acylated glucosamine unit (simplified)
    n_acyl_glucosamine_smarts = "[C@@H]1([C@@H](O)[C@H](O)[C@H](O[C@H]1CO)NC(=O)[C,C])"
    n_acyl_glucosamine = Chem.MolFromSmarts(n_acyl_glucosamine_smarts)
    if not mol.HasSubstructMatch(n_acyl_glucosamine):
        return False, "No N-acylated glucosamine units found"

    # Define pattern for KDO unit (simplified)
    kdo_smarts = "O=C([O-])[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O"
    kdo = Chem.MolFromSmarts(kdo_smarts)
    if not mol.HasSubstructMatch(kdo):
        return False, "No KDO (octulosonic acid) units found"

    # Define pattern for heptose unit (simplified)
    heptose_smarts = "O[C@H]1[C@@H](O)[C@H](O)[C@H](O)[C@H](O)[C@H](CO)O1"
    heptose = Chem.MolFromSmarts(heptose_smarts)
    heptose_matches = mol.GetSubstructMatches(heptose)
    if len(heptose_matches) < 2:
        return False, f"Found {len(heptose_matches)} heptose units, need at least 2"

    # Check for sugar rings (oxygen-containing rings)
    ring_info = mol.GetRingInfo()
    num_rings = len(ring_info.AtomRings())
    if num_rings < 5:
        return False, f"Too few rings ({num_rings}) to be a lipopolysaccharide"

    # Check for presence of long-chain fatty acids (e.g., tetradecanoic acid chains)
    # Simplified pattern for fatty acid chain (C12 or longer)
    fatty_acid_pattern = Chem.MolFromSmarts("C(=O)[C@H](O)CCCCCCCCCCCC")  # C12 chain
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 1:
        return False, "No long-chain fatty acids found"

    return True, "Molecule has N-acylated glucosamine units, KDO, heptoses, and sufficient sugar rings indicative of a lipopolysaccharide"