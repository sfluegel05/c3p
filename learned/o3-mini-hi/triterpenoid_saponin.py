"""
Classifies: CHEBI:61778 triterpenoid saponin
"""
"""
Classifies: A terpene glycoside in which the terpene moiety is a triterpenoid,
i.e. a triterpenoid saponin.
The classification checks for the combination of a sugar (glycosidic) portion and
a triterpenoid-like aglycone.
The aglycone is expected to have roughly 25-35 carbons, 4-6 rings, and no aromatic rings.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_triterpenoid_saponin(smiles: str):
    """
    Determines if a molecule is a triterpenoid saponin.
    A triterpenoid saponin is defined as a glycoside in which the non-sugar (aglycone)
    part is a triterpenoid (usually ~30 carbons, 4-6 rings, and non-aromatic).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is a triterpenoid saponin, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # First, check for evidence of at least one sugar ring.
    # We use two approximate SMARTS patterns for pyranose and furanose rings.
    # (These patterns are not 100% general but should capture many sugars.)
    sugar_smarts_list = [
        "[C@@H]1OC(CO)C(O)C(O)C1O",  # common hexopyranose (like glucose)
        "[C@@H]1OC(CO)C(O)C1O"         # a slightly smaller ring (furanose-like)
    ]
    sugar_matches = []
    for smarts in sugar_smarts_list:
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is not None:
            matches = mol.GetSubstructMatches(pattern)
            if matches:
                sugar_matches.extend(matches)
    if not sugar_matches:
        return False, "No typical sugar moiety found in molecule."
    
    # Now try to remove all sugar-like substructures.
    # We iterate through our list of sugar SMARTS and delete those substructures.
    # Note: Chem.DeleteSubstructs returns a new molecule with all the matches removed.
    aglycone_candidate = Chem.Mol(mol)  # copy original molecule
    for smarts in sugar_smarts_list:
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is not None:
            # Delete all substructures matching the sugar pattern.
            aglycone_candidate = Chem.DeleteSubstructs(aglycone_candidate, pattern)
    # Clean up: sometimes deletion leaves extra dummy atoms (e.g. "*")
    aglycone_candidate = Chem.RemoveHs(aglycone_candidate)
    # Get individual fragments from the remaining molecule.
    fragments = Chem.GetMolFrags(aglycone_candidate, asMols=True, sanitizeFrags=True)
    if not fragments:
        return False, "No aglycone fragment found after sugar removal."
    
    # Select the largest fragment by carbon count as the aglycone candidate.
    def carbon_count(m):
        return sum(1 for atom in m.GetAtoms() if atom.GetAtomicNum() == 6)
    aglycone = max(fragments, key=carbon_count)
    aglycone_carbons = carbon_count(aglycone)
    
    # Calculate number of rings (using SSSR)
    ring_count = Chem.rdMolDescriptors.CalcNumRings(aglycone)
    
    # Check for aromatic atoms; triterpenoids are aliphatic
    if any(atom.GetIsAromatic() for atom in aglycone.GetAtoms()):
        return False, "Aglycone contains aromatic rings, not typical of triterpenoids."
    
    # Define acceptable ranges (these thresholds are approximate)
    if not (25 <= aglycone_carbons <= 35):
        return False, f"Aglycone carbon count ({aglycone_carbons}) is not in the expected range (25-35)."
    if not (4 <= ring_count <= 6):
        return False, f"Aglycone ring count ({ring_count}) is not in the expected range (4-6)."
    
    # If all criteria are met, return positive classification.
    reason = (f"Molecule contains a sugar moiety and a triterpenoid-like aglycone "
              f"(aglycone carbons: {aglycone_carbons}, rings: {ring_count}).")
    return True, reason

# Example usage (uncomment to test):
# result, explanation = is_triterpenoid_saponin("OC1(C2C=3...etc")
# print(result, explanation)