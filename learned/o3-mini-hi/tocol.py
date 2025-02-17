"""
Classifies: CHEBI:39437 tocol
"""
"""
Classifies: tocol
Definition: A chromanol (2,3-dihydro-1-benzopyran) with a chroman-6-ol skeleton that is substituted at position 2 
            by a saturated or triply-unsaturated hydrocarbon chain (three isoprenoid units).
            
This implementation uses rdkit to:
  1. Find a chroman core (a 2,3-dihydro-1-benzopyran) using a SMARTS pattern.
  2. Verify that the aromatic portion of the core has a –OH substituent.
  3. Identify the substituent at the 2‑position (the carbon immediately following the heterocyclic oxygen in the saturated ring).
  4. Traverse and extract the side chain from that point.
  5. Check that the side chain is aliphatic (only carbons), has a minimum number of carbons (heuristically ≥ 10),
     and that it is either saturated (0 double bonds) or has exactly 3 double bonds.
     
If any of the checks fail, the function returns False along with a textual explanation.
Note: Tocol chemistry is complex; this heuristic may miss or misclassify some edge cases.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.rdchem import BondType

def is_tocol(smiles: str):
    """
    Determines whether a molecule given as a SMILES string belongs to the tocol class
    (a chromanol core with a specific 2‑position side chain as defined).

    Args:
       smiles (str): SMILES string of the molecule.
       
    Returns:
       bool: True if the molecule meets the tocol definition, False otherwise.
       str: Explanation for the decision.
    """
    # Parse the SMILES string into an RDKit Mol object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Step 1. Identify the chromanol (benzopyran) core.
    # We use a SMARTS that matches a 2,3-dihydro-1-benzopyran scaffold.
    # The pattern "O1CCc2ccccc2C1" matches a saturated oxygen-containing six-membered ring
    # fused with a benzene ring.
    core_smarts = "O1CCc2ccccc2C1"
    core_query = Chem.MolFromSmarts(core_smarts)
    core_matches = mol.GetSubstructMatches(core_query, useChirality=False)
    if not core_matches:
        return False, "Chromanol core (benzopyran) not found"
    
    # For our purposes we select the first match.
    core_match = core_matches[0]
    core_atom_indices = set(core_match)
    
    # Step 2. Check for the presence of a phenolic -OH on the aromatic (benzene) portion.
    # We check for an aromatic carbon directly bonded to an -OH.
    # (This is a heuristic for the chroman-6-ol part.)
    phenol_query = Chem.MolFromSmarts("c[OH]")
    if not mol.HasSubstructMatch(phenol_query):
        return False, "Aromatic hydroxyl group (phenolic OH) not found on the core"
    
    # Step 3. Locate the substituent at the 2‑position.
    # By convention in the SMARTS "O1CCc2ccccc2C1",
    # the saturated ring atoms are ordered [O, C, C, ...].
    # We take the carbon at index 1 (the first C after the O) as the position 2.
    pos2_idx = core_match[1]
    pos2_atom = mol.GetAtomWithIdx(pos2_idx)
    
    # The side chain should be attached at position 2. We find a neighbor of pos2_atom that is NOT in the core.
    side_chain_root = None
    for nbr in pos2_atom.GetNeighbors():
        if nbr.GetIdx() not in core_atom_indices:
            side_chain_root = nbr
            break
    if side_chain_root is None:
        return False, "Side chain at the 2‑position not found"
    
    # Step 4. Extract the side chain by traversing outwards from side_chain_root,
    # avoiding any atoms that are part of the core.
    side_chain_idxs = set()
    to_visit = [side_chain_root.GetIdx()]
    while to_visit:
        current_idx = to_visit.pop()
        if current_idx in side_chain_idxs:
            continue
        side_chain_idxs.add(current_idx)
        current_atom = mol.GetAtomWithIdx(current_idx)
        for nbr in current_atom.GetNeighbors():
            if nbr.GetIdx() in core_atom_indices:
                continue
            if nbr.GetIdx() not in side_chain_idxs:
                to_visit.append(nbr.GetIdx())
    
    if not side_chain_idxs:
        return False, "Failed to extract the side chain at the 2‑position"
    
    # Step 5. Verify that the side chain is aliphatic (contains only carbon atoms).
    for idx in side_chain_idxs:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() != 6:
            return False, f"Side chain contains a non-carbon atom: {atom.GetSymbol()} (atom idx: {idx})"
    
    n_chain_carbons = len(side_chain_idxs)
    if n_chain_carbons < 10:
        return False, f"Side chain is too short (found {n_chain_carbons} carbons; expected at least 10)"
    
    # Step 6. Count the number of double bonds in the side chain.
    double_bond_count = 0
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtomIdx()
        a2 = bond.GetEndAtomIdx()
        # Count only bonds where both atoms belong to the side chain.
        if a1 in side_chain_idxs and a2 in side_chain_idxs:
            if bond.GetBondType() == BondType.DOUBLE:
                double_bond_count += 1
    
    # Accept either a fully saturated side chain (0 double bonds) or one with exactly three double bonds.
    if double_bond_count not in (0, 3):
        return False, f"Side chain has {double_bond_count} double bond(s); expected 0 (fully saturated) or 3 (triply unsaturated)"
    
    return True, "Molecule contains a chromanol core with a phenolic -OH and a 2‑position side chain matching the tocol definition"

# For testing (remove or modify these lines as needed):
if __name__ == "__main__":
    # Test with a few example SMILES from the provided list.
    test_smiles_list = [
        "CC(C)CCC[C@H](C)CCC[C@@H](C)CCC[C@@]1(C)CCC2=C(C)C(O)=C(C)C(C)=C2O1",  # (S,R,S)-alpha-tocopherol
        "CC1=C(C(=C2CCC(OC2=C1C)(C)CCCC(C)CCCC(C)CCCC(C)C)C)O",  # 2,5,7,8-tetramethyl-... (tocochroman)
        "CC(C)=CCC\\C(C)=C\\CC\\C(C)=C\\CC[C@]1(C)CCc2c(C)c(O)c(C)c(C)c2O1",  # alpha-tocotrienol 
    ]
    for s in test_smiles_list:
        flag, reason = is_tocol(s)
        print("SMILES:", s)
        print("Result:", flag, "|", reason)
        print("------------------------------------------------")