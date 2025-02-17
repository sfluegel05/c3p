"""
Classifies: CHEBI:59549 essential fatty acid
"""
"""
Classifies: Essential fatty acid 
Definition: A free (non-esterified) fatty acid that is acyclic, has exactly one terminal carboxylic acid group,
a sufficiently long (≥10 carbon) chain, at least 2 carbon–carbon double bonds (apart from the acid carbonyl), 
and no extra heteroatoms beyond the atoms in a simple fatty acid (C, H, and exactly the 2 oxygens of the acid).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_essential_fatty_acid(smiles: str):
    """
    Determines if a molecule is a free essential fatty acid based on its SMILES string.
    The criteria are:
      - The SMILES string must be valid.
      - The molecule must be acyclic (no rings).
      - The molecule must contain exactly one free carboxylic acid group. In a free fatty acid,
        the acid carbon must be terminal (i.e. connected to exactly one other carbon).
      - The molecule must contain exactly two oxygen atoms (those in the acid group) and no other heteroatoms.
      - The molecule must contain at least 10 carbons.
      - Excluding the carbonyl of the acid group, the molecule must have at least 2 carbon–carbon double bonds.
      
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule meets the criteria for an essential fatty acid, False otherwise.
        str: Explanation for the decision.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure molecule is acyclic.
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings, not a simple acyclic fatty acid"
    
    # Check that the molecule does not have extra heteroatoms.
    # For our purposes, only C, H, and O are allowed.
    allowed_atomic_nums = {1, 6, 8}  # H, C, O
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, "Molecule has additional heteroatoms (e.g., N, P) not found in a simple fatty acid"
    
    # Find free carboxylic acid groups.
    # Use a SMARTS pattern that matches both protonated and ionized forms.
    # Note: The pattern matches the acid carbon as the first atom.
    acid_smarts = Chem.MolFromSmarts("[CX3](=O)[O;H1,-]")
    acid_matches = mol.GetSubstructMatches(acid_smarts)
    
    # Deduplicate acid carbons (first atom in the matched tuple).
    acid_carbon_indices = set(match[0] for match in acid_matches)
    if len(acid_carbon_indices) != 1:
        return False, f"Expected one free carboxylic acid group, found {len(acid_carbon_indices)}"
    
    # Check that the acid carbon is terminal (i.e. it has exactly one carbon neighbor).
    acid_carbon = mol.GetAtomWithIdx(next(iter(acid_carbon_indices)))
    carbon_neighbors = [nbr for nbr in acid_carbon.GetNeighbors() if nbr.GetAtomicNum() == 6]
    if len(carbon_neighbors) != 1:
        return False, "The carboxylic acid group is not terminal (acid carbon should have exactly one carbon neighbor)"
    
    # Verify total oxygen count: for a simple free fatty acid we expect exactly 2 (both in the acid group).
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count != 2:
        return False, f"Molecule contains {oxygen_count} oxygen atoms; expected exactly 2 for a free fatty acid"
    
    # Count the total number of carbon atoms.
    carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    num_carbons = len(carbons)
    if num_carbons < 10:
        return False, f"Not enough carbon atoms (found {num_carbons}, need at least 10)"
    
    # Count the number of carbon–carbon double bonds that are not part of the acid carbonyl.
    # We ignore any double bond that involves the acid carbon.
    cc_double_bonds = 0
    acid_carbon_idx = acid_carbon.GetIdx()
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
                if a1.GetIdx() == acid_carbon_idx or a2.GetIdx() == acid_carbon_idx:
                    continue
                cc_double_bonds += 1
    if cc_double_bonds < 2:
        return False, f"Not polyunsaturated enough (found {cc_double_bonds} carbon–carbon double bond(s); need at least 2)"
    
    return True, (f"Contains a terminal free carboxylic acid group, {num_carbons} carbons, and {cc_double_bonds} carbon–carbon double bonds, "
                  "consistent with an essential fatty acid")