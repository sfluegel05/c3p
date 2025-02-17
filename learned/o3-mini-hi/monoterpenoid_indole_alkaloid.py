"""
Classifies: CHEBI:65323 monoterpenoid indole alkaloid
"""
"""
Classifies: monoterpenoid indole alkaloid
A monoterpenoid indole alkaloid is biosynthesised from L-tryptophan and a diisoprenoid (usually secologanin)
building block. The heuristic approach:
  1. Look for an indole core:
       • First try a free (or N-substituted) indole via SMARTS.
       • If not found, look for a fused ring system where a 5-membered ring containing at least one nitrogen
         shares ≥2 atoms with a fully aromatic 6-membered ring.
  2. Look for a terpene-like fragment. Instead of any non‐aromatic vinyl bond, we require that at least one non‐aromatic 
     C=C bond (vinyl fragment) is “attached” to the indole by being near (within 4 bonds) to any atom of the indole core.
  3. Basic size check (minimum number of carbons and molecular weight).
  
Note: This heuristic will not capture every nuance of biosynthesis.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monoterpenoid_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid indole alkaloid using heuristic substructure patterns.
    
    The classification is based on:
      - The presence of an indole core (free, substituted, or detected via ring fusion).
      - The presence of a non‐aromatic vinyl (C=C) bond that is found in close proximity (within 4 bonds)
        to the indole core (serving as a proxy for a terpene attachment).
      - Basic size requirements in terms of carbon count and molecular weight.
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a monoterpenoid indole alkaloid, False otherwise.
        str: A reason describing the classification outcome.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # ----- 1. Identify the indole core -----
    # Try SMARTS for free indole and N-substituted indole.
    indole_free = Chem.MolFromSmarts("c1ccc2c(c1)[nH]c2")
    indole_sub  = Chem.MolFromSmarts("c1ccc2c(c1)[n]c2")
    
    indole_atoms = None
    if mol.HasSubstructMatch(indole_free):
        indole_atoms = set(mol.GetSubstructMatch(indole_free))
    elif mol.HasSubstructMatch(indole_sub):
        indole_atoms = set(mol.GetSubstructMatch(indole_sub))
    else:
        # If simple SMARTS did not match, try to find a fused ring system.
        rings = mol.GetRingInfo().AtomRings()
        for ring1 in rings:
            if len(ring1) == 5:
                # Check if the 5-membered ring has at least one nitrogen
                if any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 7 for idx in ring1):
                    for ring2 in rings:
                        if len(ring2) == 6:
                            # Check that every atom in the 6-membered ring is aromatic.
                            if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring2):
                                # Check for fusion: share at least 2 atoms.
                                if len(set(ring1).intersection(ring2)) >= 2:
                                    indole_atoms = set(ring1).union(ring2)
                                    break
                    if indole_atoms is not None:
                        break
    if indole_atoms is None:
        return False, "No indole core (free, substituted, or fused) found"
    
    # ----- 2. Look for terpene-like (isoprenoid) fragment -----
    # We iterate over all bonds; for each non-aromatic C=C bond, we check if it is "attached" to the indole core.
    terpene_found = False
    dm = Chem.GetDistanceMatrix(mol)  # distance matrix for path lengths
    for bond in mol.GetBonds():
        # Check if bond is a double bond and non-aromatic.
        if bond.GetBondTypeAsDouble() == 2 and not bond.GetIsAromatic():
            i = bond.GetBeginAtomIdx()
            j = bond.GetEndAtomIdx()
            atom_i = mol.GetAtomWithIdx(i)
            atom_j = mol.GetAtomWithIdx(j)
            # We require that both atoms are carbons (typical for a vinyl fragment).
            if atom_i.GetAtomicNum() != 6 or atom_j.GetAtomicNum() != 6:
                continue
            # Now check distance to the indole core:
            # If either atom of the double bond is within 4 bonds of any atom in the indole core,
            # we consider the vinyl (terpene) fragment as attached.
            for idx in indole_atoms:
                if dm[i][idx] <= 4 or dm[j][idx] <= 4:
                    terpene_found = True
                    break
            if terpene_found:
                break
    if not terpene_found:
        return False, "No attached non‐aromatic vinyl fragment (terpenoid moiety proxy) detected"
    
    # ----- 3. Basic size and composition checks -----
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 15:
        return False, "Not enough carbon atoms to harbor both tryptophan and a monoterpenoid unit"
    
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250:
        return False, "Molecular weight too low for a typical monoterpenoid indole alkaloid"
    
    return True, "Contains an indole core and an attached terpene-related vinyl fragment with adequate size"

# (Optional) You can run tests by calling the function with example SMILES.
# e.g., print(is_monoterpenoid_indole_alkaloid("C/C=C\\1/CN2[C@]3(C[C@@]1(C(CO)(C(=O)OC)[C@@]45C[C@@]2(O[C@@]34NC6=CC=CC=C65)[H])[H])[H]"))