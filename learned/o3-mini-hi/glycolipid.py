"""
Classifies: CHEBI:33563 glycolipid
"""
"""
Classifies: glycolipid
Definition: Any member of the class of 1,2-di-O-acylglycerols joined at oxygen 3 by a glycosidic linkage
to a carbohydrate part (usually a mono-, di- or tri-saccharide). Some bacterial glycolipids have the sugar part acylated,
and the glycerol component may be absent.
This program uses several heuristic tests:
  1. Identify carbohydrate rings as non‐aromatic 5- or 6‐membered rings containing exactly one oxygen.
  2. Confirm that at least one sugar ring is linked via an ether (glycosidic) bond (from one of its carbons) to an external fragment.
  3. Require the presence of an acyl linkage (ester or amide) that could anchor a fatty acid.
  4. Check that there is a long carbon chain—now measured over all non‐aromatic, non‐ring carbon atoms regardless of hybridization—of at least 8 atoms.
  5. Ensure that the molecular weight is high enough (>=400 Da) and that sufficient rotatable bonds (>=3) are present.
  6. Filter out molecules containing an adenine substructure (indicative of CoA/nucleotide cofactors) that are not glycolipids.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import AllChem

def get_sugar_rings(mol):
    """
    Identify sugar rings: non-aromatic ring of size 5 or 6 having exactly one oxygen.
    Returns a list of sets (each set is atom indices in that ring).
    """
    sugar_rings = []
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) not in [5,6]:
            continue
        # Exclude aromatic rings.
        if any(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            continue
        oxy_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
        if oxy_count == 1:
            sugar_rings.append(set(ring))
    return sugar_rings

def has_glycosidic_linkage(mol, sugar_rings):
    """
    Check if at least one sugar ring is connected via an ether bond (through an oxygen)
    from one of its carbons to a non-sugar fragment.
    """
    for ring in sugar_rings:
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # We look at carbons in the ring.
            if atom.GetAtomicNum() != 6:
                continue
            for nbr in atom.GetNeighbors():
                # Looking for an oxygen not in the sugar ring.
                if nbr.GetAtomicNum() == 8 and (nbr.GetIdx() not in ring):
                    # Check that the oxygen is further connected to a non-sugar atom.
                    for nn in nbr.GetNeighbors():
                        if nn.GetIdx() != atom.GetIdx() and (nn.GetIdx() not in ring):
                            return True
    return False

def longest_nonring_carbon_chain(mol):
    """
    Computes the length of the longest chain of connected carbon atoms that are:
      - Carbon atoms (atomic number 6)
      - Not aromatic
      - Not in any ring.
    Unlike our previous version we no longer restrict to sp3 centers so unsaturated (alkene) chains count.
    We use a DFS search over eligible atoms.
    """
    eligible = set()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and (not atom.GetIsAromatic()) and (not atom.IsInRing()):
            eligible.add(atom.GetIdx())
    neighbor_map = {}
    for idx in eligible:
        atom = mol.GetAtomWithIdx(idx)
        neighbor_map[idx] = [nbr.GetIdx() for nbr in atom.GetNeighbors() if nbr.GetIdx() in eligible]
    max_chain = 0
    def dfs(node, visited):
        current_length = 1
        for nbr in neighbor_map.get(node, []):
            if nbr not in visited:
                current_length = max(current_length, 1 + dfs(nbr, visited | {nbr}))
        return current_length
    for idx in eligible:
        chain_length = dfs(idx, {idx})
        if chain_length > max_chain:
            max_chain = chain_length
        if max_chain >= 8:
            break
    return max_chain

def is_glycolipid(smiles: str):
    """
    Determines if a molecule is a glycolipid based on its SMILES string.
    Heuristic tests:
      1. Identify at least one non-aromatic sugar ring (5- or 6-membered with one oxygen).
      2. Verify that a sugar ring is connected via an ether (glycosidic) bond to an external fragment.
      3. Detect acyl linkages by looking for ester ([OX2][CX3](=O)[#6]) or amide (N[C;!R](=O)[#6]) patterns.
      4. Ensure there is a long acyl or non-ring carbon chain (>=8 carbons; may be saturated or unsaturated).
      5. Check overall molecular weight (>=400 Da) and at least 3 rotatable bonds.
      6. Filter out molecules with adenine substructures, as these are typical of nucleotide or CoA derivatives.
    Returns:
        (bool, str): A tuple where the first element indicates if the molecule is classified as a glycolipid,
                     and the second gives the reason for its classification or why it failed.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # (A) Filter out molecules containing an adenine moiety (common in CoA/nucleotide derivatives).
    adenine_smarts = Chem.MolFromSmarts("c1ncnc2ncnc(c12)")
    if adenine_smarts and mol.HasSubstructMatch(adenine_smarts):
        return False, "Molecule contains an adenine moiety (likely a nucleotide/coenzyme derivative), not a glycolipid"

    # (1) Detect sugar rings.
    sugar_rings = get_sugar_rings(mol)
    if not sugar_rings:
        return False, "No carbohydrate (sugar ring) moiety found"

    # (2) Check for a glycosidic linkage.
    if not has_glycosidic_linkage(mol, sugar_rings):
        return False, "No glycosidic linkage (ether bond between a sugar ring and non-sugar fragment) found"

    # (3) Look for acyl linkages that could anchor fatty acid chains.
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=O)[#6]")
    amide_pattern = Chem.MolFromSmarts("N[C;!R](=O)[#6]")
    ester_matches = mol.GetSubstructMatches(ester_pattern) if ester_pattern is not None else []
    amide_matches = mol.GetSubstructMatches(amide_pattern) if amide_pattern is not None else []
    if (len(ester_matches) + len(amide_matches)) < 1:
        return False, "No acyl linkage (ester or amide) found that could anchor a fatty acyl chain"

    # (4) Check for a long aliphatic (or unsaturated) carbon chain.
    chain_len = longest_nonring_carbon_chain(mol)
    if chain_len < 8:
        return False, f"No long carbon chain detected (longest chain has only {chain_len} carbons; need >= 8)"

    # (5) Check molecular weight and rotatable bonds.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for a glycolipid"
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 3:
        return False, "Not enough rotatable bonds to support a long fatty acyl chain"

    return True, "Molecule contains a sugar ring linked via a glycosidic bond and an acyl link to a long carbon chain, consistent with a glycolipid structure"

# For testing, you may uncomment the following lines:
# if __name__ == "__main__":
#     test_smiles = "OC[C@H]1O[C@@H](O[C@@H]2[C@@H](O)[C@H](O)[C@@H](CO)O[C@H]2OCCCCCCCC\\C=C/CCCCCCCC(O)=O"
#     result, reason = is_glycolipid(test_smiles)
#     print("Is glycolipid?", result)
#     print("Reason:", reason)