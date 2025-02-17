"""
Classifies: CHEBI:65323 monoterpenoid indole alkaloid
"""
"""
Classifies: monoterpenoid indole alkaloid
A terpenoid indole alkaloid is biosynthesised from L-tryptophan and a diisoprenoid (usually secologanin) building blocks.
The classification below uses a heuristic approach:
  1. Check for an indole core (free or N-substituted).
  2. Look for a non‐aromatic vinyl (C=C) fragment outside of the indole ring as a proxy for the terpenoid moiety.
  3. Verify that the molecule is of sufficient size (a minimum number of carbons and overall molecular weight).
  
Note: This heuristic may not capture every nuance of biosynthesis.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monoterpenoid_indole_alkaloid(smiles: str):
    """
    Determines if the molecule is a monoterpenoid indole alkaloid using heuristic substructure patterns.
    
    The approach is based on the idea that such molecules must contain:
      - An indole core (which can be either unsubstituted at nitrogen or N-methylated/substituted).
      - A terpene-related portion. Here we check for a non‐aromatic vinyl (C=C) fragment that is not part
        of an aromatic system.
      - A minimum carbon count and molecular weight (to roughly ensure the presence of both tryptophan and
        a monoterpenoid unit).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a monoterpenoid indole alkaloid, False otherwise.
        str: Reason for the classification.
    """
    
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # --- 1. Check for the indole core ---
    # Define two SMARTS patterns: one for a free indole and one for an N-substituted indole.
    indole_smarts_free = "c1ccc2c(c1)[nH]c(c2)"       # indole with a free NH
    indole_smarts_sub = "c1ccc2c(c1)nc(c2)"             # indole with substituted nitrogen
    indole_free = Chem.MolFromSmarts(indole_smarts_free)
    indole_sub  = Chem.MolFromSmarts(indole_smarts_sub)
    
    if not (mol.HasSubstructMatch(indole_free) or mol.HasSubstructMatch(indole_sub)):
        return False, "No indole core (free or substituted) found"
        
    # --- 2. Look for terpene-like (isoprenoid) fragment ---
    # Here we use the presence of at least one non‐aromatic vinyl bond ([CH2]=[CH]) as a proxy.
    vinyl_smarts = "[CH2]=[CH]"
    vinyl = Chem.MolFromSmarts(vinyl_smarts)
    vinyl_matches = mol.GetSubstructMatches(vinyl)
    terpene_found = False
    for match in vinyl_matches:
        # Check that at least one atom in the match is not aromatic,
        # so that the double bond is unlikely to be part of an aromatic ring.
        if any(not mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in match):
            terpene_found = True
            break
    
    if not terpene_found:
        return False, "No non‐aromatic vinyl fragment (terpenoid moiety proxy) detected"
    
    # --- 3. Basic size and composition checks ---
    # Count carbon atoms.
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 15:
        return False, "Not enough carbon atoms to harbor both tryptophan and a monoterpenoid segment"
    
    # Check molecular weight (rough heuristic - many monoterpenoid indole alkaloids are >250 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250:
        return False, "Molecular weight too low for a typical monoterpenoid indole alkaloid"
    
    # If all tests pass
    return True, "Contains an indole core and a terpene-related vinyl fragment with adequate size"