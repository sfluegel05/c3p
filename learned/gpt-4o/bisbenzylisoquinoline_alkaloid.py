"""
Classifies: CHEBI:133004 bisbenzylisoquinoline alkaloid
"""
from rdkit import Chem
from rdkit.Chem import rdqueries

def is_bisbenzylisoquinoline_alkaloid(smiles: str):
    """
    Determine if a molecule is a bisbenzylisoquinoline alkaloid based on its SMILES string.
   
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule matches the bisbenzylisoquinoline alkaloid pattern, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Smart and flexible pattern for benzylisoquinoline core, accounts for possible variations
    benzylisoquinoline_pattern = Chem.MolFromSmarts('c1ccc2c(c1)CN3C=CC=CC3=C2')  # More common structure
    matches = mol.GetSubstructMatches(benzylisoquinoline_pattern)
    if len(matches) < 2:
        return False, "Does not contain two benzylisoquinoline units"

    # More detailed ether bridge checking, ensuring the link between two units
    # Look for an 'O' linking two aromatic cycles representing benzylisoquinoline structures
    ether_pattern = Chem.MolFromSmarts('c-c-O-c-c')
    ether_matches = mol.GetSubstructMatches(ether_pattern)
    if len(ether_matches) == 0:
        return False, "No specific ether bridge detected between benzylisoquinoline units"

    # Optional: look for methylenedioxy bridging between the units
    methylenedioxy_pattern = Chem.MolFromSmarts('COC1COC1')
    methylenedioxy_matches = mol.GetSubstructMatches(methylenedioxy_pattern)
    if len(methylenedioxy_matches) > 0:
        return True, "Contains methylenedioxy bridge indicating bisbenzylisoquinoline alkaloid"

    # Optional: for further confidence, additional carbon-carbon direct bridging
    # Using a more relaxed search for possible C-C bridging between benzylisoquinoline units
    cc_bridge_pattern = Chem.MolFromSmarts('c1cc[c,C]ccc1-c2cc[c,C]ccc2')
    cc_bridge_matches = mol.GetSubstructMatches(cc_bridge_pattern)
    if len(cc_bridge_matches) > 0:
        return True, "Contains C-C bridge between benzylisoquinoline units"

    return True, "Classified as bisbenzylisoquinoline alkaloid based on structural patterns"