"""
Classifies: CHEBI:33563 glycolipid
"""
"""
Classifies: glycolipid
Definition: Any member of the class of 1,2-di-O-acylglycerols joined at oxygen 3 by a glycosidic linkage
to a carbohydrate part (usually a mono-, di- or tri-saccharide).
Note: Some bacterial glycolipids have the sugar part acylated and the glycerol part may be absent.
This program uses an improved series of heuristic tests:
  1. Identify carbohydrate rings as non‐aromatic 5- or 6‐membered rings containing exactly one oxygen.
  2. Confirm that at least one sugar ring is linked via an ether (glycosidic) bond (from one of its carbons) to a non‐sugar fragment.
  3. Detect acyl linkages by looking for ester ([OX2][CX3](=O)[#6]) or amide (N[C;!R](=O)[#6]) patterns.
     Also, allow for the presence of a polyisoprenyl chain (repeating isoprene units) as an alternative lipid anchor.
  4. Require the presence of a long non‐ring carbon chain (using a SMARTS pattern requiring at least 8 contiguous aliphatic carbons,
     thus capturing both saturated and unsaturated fatty acyl chains).
  5. Check that the molecule’s overall weight (>=400 Da) and rotatable bonds (>=3) are adequate.
  6. Filter out molecules containing adenine (typical of nucleotide/coenzyme derivatives).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import AllChem

def get_sugar_rings(mol):
    """
    Identify sugar rings as non‐aromatic rings of size 5 or 6 with exactly one oxygen.
    Returns a list of sets (each set contains atom indices in that ring).
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
    Check that at least one sugar ring is connected via an ether (C–O) bond 
    (from a carbon in the ring to an oxygen outside the ring).
    """
    for ring in sugar_rings:
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6:
                continue  # Only consider carbons in the ring.
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 8 and (nbr.GetIdx() not in ring):
                    # Verify that the oxygen attaches to at least one atom outside the ring.
                    for nn in nbr.GetNeighbors():
                        if nn.GetIdx() != atom.GetIdx() and (nn.GetIdx() not in ring):
                            return True
    return False

def has_acyl_or_polyprenyl(mol):
    """
    Return True if the molecule contains either an acyl (ester or amide) linkage or a polyisoprenyl (repeating isoprene)
    motif.
    """
    # Look for ester and amide bonds:
    ester_smarts = "[OX2][CX3](=O)[#6]"
    amide_smarts = "N[C;!R](=O)[#6]"
    ester_pattern = Chem.MolFromSmarts(ester_smarts)
    amide_pattern = Chem.MolFromSmarts(amide_smarts)
    acyl_count = 0
    if ester_pattern:
        acyl_count += len(mol.GetSubstructMatches(ester_pattern))
    if amide_pattern:
        acyl_count += len(mol.GetSubstructMatches(amide_pattern))
    # Allow an alternative if a polyisoprenyl chain is present.
    # This pattern looks for at least 2 repeating isoprene units.
    polyisoprenyl_smarts = "([CH2]C(=C)[CH2]){2,}"
    poly_pattern = Chem.MolFromSmarts(polyisoprenyl_smarts)
    has_poly = False
    if poly_pattern and mol.HasSubstructMatch(poly_pattern):
        has_poly = True

    return (acyl_count >= 1) or has_poly

def has_long_nonring_chain(mol):
    """
    Check for the existence of at least one linear (non-ring) aliphatic chain.
    We use a SMARTS query that looks for at least 8 contiguous carbon atoms that are non-aromatic and not in rings.
    """
    # This SMARTS pattern matches 8 or more consecutive aliphatic (non-ring, non-aromatic) carbon atoms.
    chain_smarts = "[C;!R;!a]{8,}"
    chain_pattern = Chem.MolFromSmarts(chain_smarts)
    if chain_pattern and mol.HasSubstructMatch(chain_pattern):
        return True
    return False

def is_glycolipid(smiles: str):
    """
    Determines if a molecule is a glycolipid based on its SMILES string.
    The classifier uses several heuristics:
      1. Identifies at least one non-aromatic sugar ring (5- or 6-membered with one oxygen).
      2. Verifies that at least one sugar ring is connected by an ether linkage to a non-sugar fragment.
      3. Detects acyl linkages (ester/amide) or a polyisoprenyl chain as a lipid anchor.
      4. Requires the presence of a long linear aliphatic chain (>=8 contiguous carbons, non‐ring).
      5. Checks that the molecular weight is >=400 Da and there are at least 3 rotatable bonds.
      6. Filters out molecules containing adenine (which indicate nucleotide/coenzyme).
      
    Returns:
        (bool, str): True and a message if classified as a glycolipid; otherwise False and a reason.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # (A) Filter out nucleotide or coenzyme derivatives by rejecting adenine motifs.
    adenine_smarts = "c1ncnc2ncnc(c12)"
    adenine_pattern = Chem.MolFromSmarts(adenine_smarts)
    if adenine_pattern and mol.HasSubstructMatch(adenine_pattern):
        return False, "Molecule contains an adenine moiety (likely a nucleotide/coenzyme derivative), not a glycolipid"

    # (1) Identify sugar rings.
    sugar_rings = get_sugar_rings(mol)
    if not sugar_rings:
        return False, "No carbohydrate (sugar ring) moiety found"

    # (2) Check for a glycosidic (ether) linkage from a sugar ring.
    if not has_glycosidic_linkage(mol, sugar_rings):
        return False, "No glycosidic linkage (ether bond between a sugar ring and non-sugar fragment) found"

    # (3) Check for an acyl linkage (ester or amide) or for a polyisoprenyl lipid anchor.
    if not has_acyl_or_polyprenyl(mol):
        return False, "No acyl linkage (ester/amide) or polyisoprenyl motif found that could anchor a fatty acyl chain"

    # (4) Check for a long fatty acyl or lipophilic chain.
    if not has_long_nonring_chain(mol):
        return False, "No long aliphatic chain detected (need at least 8 contiguous non-ring carbons)"

    # (5) Check molecular weight and rotatable bonds.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for a glycolipid"
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 3:
        return False, "Not enough rotatable bonds to support a flexible fatty acyl chain"

    return True, "Molecule contains a sugar ring linked via a glycosidic bond and a fatty acyl (or polyisoprenyl) chain, consistent with glycolipid structure"

# Example testing (uncomment for testing):
# if __name__ == "__main__":
#     test_smiles = "OC[C@H]1O[C@@H](O[C@@H]2[C@@H](O)[C@H](O)[C@@H](CO)O[C@H]2OCCCCCCCC\\C=C/CCCCCCCC(O)=O"
#     result, reason = is_glycolipid(test_smiles)
#     print("Is glycolipid?", result)
#     print("Reason:", reason)