"""
Classifies: CHEBI:65323 monoterpenoid indole alkaloid
"""
"""
Classifies: monoterpenoid indole alkaloid
A terpenoid indole alkaloid is biosynthesised from L-tryptophan and a diisoprenoid (usually secologanin) building block.
Heuristic approach:
  1. Check for an indole core. First test via standard SMARTS (for free or substituted indoles).
     If no match, search for a fused ring system where a 5-membered ring containing at least one nitrogen
     shares ≥2 atoms with a fully aromatic 6-membered ring.
  2. Look for a terpene moiety proxy by searching for any non‐aromatic vinyl (C=C) fragment,
     which is broadened to catch substituted alkenes.
  3. Check that the molecule is of sufficient size.
  
Note: This heuristic will not capture every nuance of biosynthesis.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monoterpenoid_indole_alkaloid(smiles: str):
    """
    Determines if the molecule is a monoterpenoid indole alkaloid using heuristic substructure patterns.
    
    The classification is based on:
      - The presence of an indole core (either free, substituted, or detected via fused rings).
      - A terpene-related fragment detected as a non‐aromatic vinyl (C=C) bond.
      - Basic size requirements in terms of carbon count and molecular weight.
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a monoterpenoid indole alkaloid, False otherwise.
        str: A reason describing the classification.
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # --- 1. Check for the indole core ---
    # Define SMARTS patterns for a free indole and an N-substituted indole.
    indole_free = Chem.MolFromSmarts("c1ccc2c(c1)[nH]c2")
    indole_sub  = Chem.MolFromSmarts("c1ccc2c(c1)[n]c2")
    
    found_indole = mol.HasSubstructMatch(indole_free) or mol.HasSubstructMatch(indole_sub)
    
    # If not found by the simple patterns, try to detect an indole by ring fusion.
    # Look for a 5-membered ring (with at least one nitrogen) and a 6-membered aromatic ring
    # that share at least two atoms.
    if not found_indole:
        rings = mol.GetRingInfo().AtomRings()
        for ring1 in rings:
            if len(ring1) == 5 and any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 7 for idx in ring1):
                for ring2 in rings:
                    if len(ring2) == 6:
                        # Check that every atom in the 6-membered ring is aromatic.
                        if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring2):
                            # Check for fusion (sharing of ≥2 atoms)
                            if len(set(ring1).intersection(ring2)) >= 2:
                                found_indole = True
                                break
                if found_indole:
                    break
    if not found_indole:
        return False, "No indole core (free, substituted, or fused) found"
    
    # --- 2. Look for terpene-like (isoprenoid) fragment ---
    # Use a generalized non‐aromatic vinyl bond SMARTS pattern.
    vinyl = Chem.MolFromSmarts("[C;!a]=[C;!a]")
    vinyl_matches = mol.GetSubstructMatches(vinyl)
    terpene_found = len(vinyl_matches) > 0
    if not terpene_found:
        return False, "No non‐aromatic vinyl fragment (terpenoid moiety proxy) detected"
    
    # --- 3. Basic size and composition checks ---
    # Count carbon atoms.
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 15:
        return False, "Not enough carbon atoms to harbor both tryptophan and a monoterpenoid unit"
    
    # Check molecular weight.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250:
        return False, "Molecular weight too low for a typical monoterpenoid indole alkaloid"
    
    return True, "Contains an indole core and a terpene-related vinyl fragment with adequate size"