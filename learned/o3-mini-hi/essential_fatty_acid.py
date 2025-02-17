"""
Classifies: CHEBI:59549 essential fatty acid
"""
"""
Classifies: Essential fatty acid
Definition: A free (non-esterified) fatty acid that is acyclic, has exactly one terminal carboxylic acid group,
contains no heteroatoms beyond C, H, and exactly 2 oxygens (of the acid), has a sufficiently long chain (≥16 carbons),
and has at least 2 carbon–carbon double bonds aside from the acid carbonyl.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_essential_fatty_acid(smiles: str):
    """
    Determines if a molecule is a free essential fatty acid based on its SMILES string.
    Criteria:
      - Valid SMILES.
      - Acyclic (no rings).
      - Contains exactly one free carboxylic acid group (carboxyl carbon is terminal, i.e. has one C neighbor).
      - Contains only C, H, and O; exactly two oxygen atoms overall.
      - Contains a sufficiently long chain: at least 16 carbon atoms.
      - Excluding the acid carbonyl, the molecule must have at least 2 carbon–carbon double bonds.
      
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule meets the criteria, False otherwise.
        str: Explanation for the decision.
    """
    # Convert the SMILES string to a molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check that the molecule is acyclic.
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings, not a simple acyclic fatty acid"
    
    # Ensure the molecule does not have extra heteroatoms. Only C, H, and O are allowed.
    allowed_atomic_nums = {1, 6, 8}  # H, C, O
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, "Molecule has additional heteroatoms (e.g., N, P) not found in a simple fatty acid"
    
    # Locate the free carboxylic acid group using SMARTS. This pattern matches both protonated and ionized forms.
    acid_smarts = Chem.MolFromSmarts("[CX3](=O)[O;H1,-]")
    acid_matches = mol.GetSubstructMatches(acid_smarts)
    # Deduplicate acid carbons (first atom in the match tuple).
    acid_carbon_indices = set(match[0] for match in acid_matches)
    if len(acid_carbon_indices) != 1:
        return False, f"Expected one free carboxylic acid group, found {len(acid_carbon_indices)}"
    
    # Check that the acid carbon is terminal.
    acid_carbon = mol.GetAtomWithIdx(next(iter(acid_carbon_indices)))
    carbon_neighbors = [nbr for nbr in acid_carbon.GetNeighbors() if nbr.GetAtomicNum() == 6]
    if len(carbon_neighbors) != 1:
        return False, "The carboxylic acid group is not terminal (acid carbon should have exactly one carbon neighbor)"
    
    # Verify the total oxygen count: for a simple free fatty acid we expect exactly 2 (from the acid group).
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count != 2:
        return False, f"Molecule contains {oxygen_count} oxygen atoms; expected exactly 2 for a free fatty acid"
    
    # Count the total number of carbon atoms.
    carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    num_carbons = len(carbons)
    if num_carbons < 16:
        return False, f"Not enough carbon atoms (found {num_carbons}, need at least 16)"
    
    # Count the number of C–C double bonds that are not part of the acid carbonyl.
    cc_double_bonds = 0
    acid_carbon_idx = acid_carbon.GetIdx()
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
                # Exclude bonds involving the acid carbon (i.e. the carbonyl)
                if a1.GetIdx() == acid_carbon_idx or a2.GetIdx() == acid_carbon_idx:
                    continue
                cc_double_bonds += 1
    if cc_double_bonds < 2:
        return False, f"Not polyunsaturated enough (found {cc_double_bonds} C–C double bond(s); need at least 2)"
    
    return True, (f"Contains a terminal free carboxylic acid group, {num_carbons} carbons, and {cc_double_bonds}"
                  " carbon–carbon double bonds, consistent with an essential fatty acid")