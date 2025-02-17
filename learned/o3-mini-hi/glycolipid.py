"""
Classifies: CHEBI:33563 glycolipid
"""
"""
Classifies: glycolipid
Definition: Any member of the class of 1,2-di-O-acylglycerols joined at oxygen 3 by a glycosidic linkage
to a carbohydrate part (usually a mono-, di- or tri-saccharide). Some bacterial glycolipids have the sugar part acylated,
and the glycerol component may be absent.
This program uses several heuristic tests:
  1. It identifies sugar rings as non-aromatic 5- or 6-membered rings containing exactly one oxygen.
  2. It confirms that at least one sugar ring is linked via an ether (glycosidic) bond (connecting a sugar ring atom to a non-sugar fragment).
  3. It requires an acyl linkage (ester or amide) that could anchor a fatty acid.
  4. It computes the longest chain of connected sp3 (non-ring) carbons to verify the presence of a long aliphatic fatty acyl chain.
  5. It also checks that the molecular weight is >500 Da and that at least 3 rotatable bonds exist.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import AllChem

def get_sugar_rings(mol):
    """
    Identify sugar rings: non-aromatic rings of size 5 or 6 having exactly one oxygen.
    Returns a list of sets (each set is atom indices in that ring).
    """
    sugar_rings = []
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) not in [5, 6]:
            continue
        # Skip aromatic rings.
        if any(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            continue
        oxy_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
        if oxy_count == 1:
            sugar_rings.append(set(ring))
    return sugar_rings

def has_glycosidic_linkage(mol, sugar_rings):
    """
    Check if at least one sugar ring is connected via an ether bond (an O atom bridging)
    from one of its (typically carbon) atoms to a non-sugar part.
    Instead of a fixed SMARTS we examine the neighbors.
    """
    for ring in sugar_rings:
        # Look at each atom in the sugar ring. For glycolipid, the linkage is usually via a carbon.
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # We are interested if this carbon is a potential anomeric carbon; here we at least require it is carbon.
            if atom.GetAtomicNum() != 6:
                continue
            # Check the neighbors of this ring carbon.
            for nbr in atom.GetNeighbors():
                # If this neighbor is an oxygen and not part of the sugar ring then consider its other neighbor.
                if nbr.GetAtomicNum() == 8 and (nbr.GetIdx() not in ring):
                    # For a glycosidic bond, the oxygen should be connected to at least one atom not in the sugar ring.
                    for nn in nbr.GetNeighbors():
                        if nn.GetIdx() != atom.GetIdx() and (nn.GetIdx() not in ring):
                            return True
    return False

def longest_aliphatic_chain(mol):
    """
    Computes the length of the longest chain of connected sp3 carbons that are:
      - carbon atoms (atomic number 6)
      - non-aromatic
      - not in any ring.
    We use a simple DFS.
    """
    # Mark exocyclic sp3 carbons (non-ring, not aromatic) as eligible.
    eligible = set()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and not atom.GetIsAromatic() and not atom.IsInRing():
            if atom.GetHybridization() == Chem.rdchem.HybridizationType.SP3:
                eligible.add(atom.GetIdx())
    # Build a neighbor map restricted to eligible atoms.
    neighbor_map = {}
    for idx in eligible:
        atom = mol.GetAtomWithIdx(idx)
        neighbor_map[idx] = [nbr.GetIdx() for nbr in atom.GetNeighbors() if nbr.GetIdx() in eligible]

    # DFS over eligible atoms.
    visited_global = set()
    max_chain = 0

    def dfs(node, visited):
        length = 1
        for nbr in neighbor_map.get(node, []):
            if nbr not in visited:
                length = max(length, 1 + dfs(nbr, visited | {nbr}))
        return length

    for idx in eligible:
        chain_length = dfs(idx, {idx})
        if chain_length > max_chain:
            max_chain = chain_length
        # Early exit if we found a sufficiently long chain
        if max_chain >= 8:
            break
    return max_chain

def is_glycolipid(smiles: str):
    """
    Determines if a molecule is a glycolipid based on its SMILES string.
    Heuristic tests:
      1. Identify a sugar moiety: scan for a non-aromatic 5- or 6-membered ring with exactly one oxygen.
      2. Ensure at least one glycosidic (ether) linkage from a sugar ring (via one of its carbons) to a non-sugar part.
      3. Identify acyl linkages by spotting ester (or amide) patterns that can anchor fatty acid chains.
      4. Verify the existence of a long aliphatic chain (longest chain of connected aliphatic sp3 carbons >= 8 atoms).
      5. Confirm that the overall molecular weight is >500 Da and that the molecule has at least 3 rotatable bonds.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): A tuple where the first element indicates if the molecule is classified as a glycolipid,
                     and the second gives the reason for its classification or why it failed.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # (1) Find sugar rings.
    sugar_rings = get_sugar_rings(mol)
    if not sugar_rings:
        return False, "No carbohydrate (sugar ring) moiety found"
    
    # (2) Check for a glycosidic linkage (ether bond bridging a sugar ring and an external fragment).
    if not has_glycosidic_linkage(mol, sugar_rings):
        return False, "No glycosidic linkage (ether bond between a sugar ring and non-sugar fragment) found"

    # (3) Look for acyl linkages.
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=O)[#6]")
    amide_pattern = Chem.MolFromSmarts("N[C;!R](=O)[#6]")
    ester_matches = mol.GetSubstructMatches(ester_pattern) if ester_pattern else []
    amide_matches = mol.GetSubstructMatches(amide_pattern) if amide_pattern else []
    if (len(ester_matches) + len(amide_matches)) < 1:
        return False, "No acyl linkage (ester or amide) found that could anchor a fatty acid chain"

    # (4) Check for a long aliphatic chain.
    chain_len = longest_aliphatic_chain(mol)
    if chain_len < 8:
        return False, f"No long aliphatic chain detected (longest chain has only {chain_len} carbon atoms; need >=8)"
    
    # (5) Check molecular weight and rotatable bonds.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for a glycolipid"
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 3:
        return False, "Not enough rotatable bonds to support long fatty acyl chains"

    return True, "Molecule contains a sugar ring linked via a glycosidic bond and acylated linkage to a long fatty acyl chain consistent with a glycolipid structure"

# For testing purposes, you may uncomment the following:
# if __name__ == "__main__":
#     test_smiles = "OC[C@H]1O[C@@H](O[C@@H]2[C@@H](O)[C@H](O)[C@@H](CO)O[C@H]2OCCCCCCCC\\C=C/CCCCCCCC(O)=O"
#     result, reason = is_glycolipid(test_smiles)
#     print("Is glycolipid?", result)
#     print("Reason:", reason)