"""
Classifies: CHEBI:76575 monoradylglycerol
"""
"""
Classifies: monoradylglycerol
Definition: Any lipid that is glycerol bearing a single acyl, alkyl or alk-1-enyl substituent at an unspecified position.
In our algorithm we search for a glycerol-like backbone – three sp3 carbons connected in a line – where each carbon bears one oxygen substituent.
Exactly one of those three oxygens is “substituted” (i.e. it has no attached hydrogen) and the connected substituent must be carbon-based,
acyclic, and contain at least a minimal number of connected carbon atoms.
"""

from rdkit import Chem

def is_monoradylglycerol(smiles: str):
    """
    Determines if a molecule is a monoradylglycerol based on its SMILES string.
    A monoradylglycerol is defined as glycerol bearing a single acyl/alkyl/alk-1-enyl substituent.
    Here we attempt to find a free glycerol (a three-carbon chain with each carbon bearing one oxygen), where exactly
    one of the oxygens is substituted (has no hydrogen) and the attached substituent comes from a carbon chain of minimal length.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is classified as a monoradylglycerol, False otherwise.
        str: A reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Helper: Given a starting atom, traverse (DFS) over connected carbons (atomic number 6)
    # that are not part of the glycerol backbone. We also require that none is in a ring.
    def count_carbon_chain(start_atom, blocked_ids):
        visited = set()
        def dfs(atom):
            if atom.GetAtomicNum() != 6:
                return 0
            # if any atom in the fragment is in a ring, we abort to avoid sugar-like fragments.
            if atom.IsInRing():
                return -100  # mark as invalid fragment
            visited.add(atom.GetIdx())
            max_count = 1
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in blocked_ids or nbr.GetIdx() in visited:
                    continue
                if nbr.GetAtomicNum() == 6:
                    cnt = 1 + dfs(nbr)
                    if cnt > max_count:
                        max_count = cnt
            return max_count
        return dfs(start_atom)

    # The minimal carbon chain length we demand for the substituent.
    # (Acetate has two carbons so we require at least 2.)
    MIN_CHAIN_LENGTH = 2

    # Iterate over atoms to find a candidate for the middle carbon of a glycerol backbone.
    # In free glycerol (or monoacyl/monoalkyl glycerol), the backbone is CH2-CHOH-CH2.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        # Only consider sp3 carbons
        if atom.GetHybridization().name != "SP3":
            continue
        # Candidate for the middle: must have exactly two sp3 carbon neighbors.
        nbr_carbons = [nbr for nbr in atom.GetNeighbors() 
                       if nbr.GetAtomicNum()==6 and nbr.GetHybridization().name=="SP3"]
        if len(nbr_carbons) != 2:
            continue
        # Also, the two neighbor carbons should not be directly bonded (to ensure a linear chain)
        if mol.GetBondBetweenAtoms(nbr_carbons[0].GetIdx(), nbr_carbons[1].GetIdx()):
            continue

        # Define backbone as [terminal1, middle, terminal2]
        # For ease of reporting, assign positions 1, 2, 3.
        backbone = [nbr_carbons[0], atom, nbr_carbons[1]]

        # Further check: For a glycerol backbone, we expect each backbone carbon to have only one oxygen neighbor (not in the backbone).
        oxygen_neighbors = []  # will be parallel to backbone so that index 0 is from position 1, etc.
        valid_backbone = True
        backbone_ids = set(a.GetIdx() for a in backbone)
        for carbon in backbone:
            # Get oxygen neighbors that are not in the backbone
            oxygens = [nbr for nbr in carbon.GetNeighbors() 
                       if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in backbone_ids]
            if len(oxygens) != 1:
                valid_backbone = False
                break
            oxygen_neighbors.append(oxygens[0])
        if not valid_backbone:
            continue

        # Count free -OH groups vs substituted oxygen.
        free_OH_count = 0
        substituted_index = -1
        for idx, o_atom in enumerate(oxygen_neighbors):
            # We rely on GetTotalNumHs: a free hydroxyl oxygen should have at least one attached H.
            # (Note: this may sometimes be implicit.)
            if o_atom.GetTotalNumHs() and o_atom.GetTotalNumHs() > 0:
                free_OH_count += 1
            else:
                substituted_index = idx
        # Exactly two free OH and one substituted oxygen required.
        if free_OH_count != 2 or substituted_index < 0:
            continue

        # Now examine the substituent attached via the substituted oxygen.
        sub_o = oxygen_neighbors[substituted_index]
        # Which backbone carbon is this oxygen attached to?
        parent_carbon = None
        for nbr in sub_o.GetNeighbors():
            if nbr.GetIdx() in backbone_ids:
                parent_carbon = nbr
                break
        if parent_carbon is None:
            continue

        # Among the oxygen's neighbors, choose the one that is not the backbone carbon.
        sub_neighbors = [nbr for nbr in sub_o.GetNeighbors() if nbr.GetIdx() not in backbone_ids]
        if not sub_neighbors:
            continue
        sub_anchor = sub_neighbors[0]
        if sub_anchor.GetAtomicNum() != 6:
            continue

        # Perform a DFS starting from sub_anchor following only carbon atoms that are NOT in the backbone.
        chain_length = count_carbon_chain(sub_anchor, blocked_ids=backbone_ids)
        # Abort if the DFS encountered a ring (indicated by negative count)
        if chain_length < MIN_CHAIN_LENGTH:
            continue

        pos = substituted_index + 1  # positions: 1 (first terminal), 2 (middle), or 3 (second terminal)

        reason = (f"Found glycerol backbone (atoms {', '.join(str(a.GetIdx()) for a in backbone)}) with substitution "
                  f"at position {pos}; substituent chain length = {chain_length}, free -OH groups = 2.")
        return True, reason

    return False, "No glycerol backbone with one substituted oxygen (and two free –OH groups) attached to a valid acyl/alkyl chain was found."


# Example usage and testing with several SMILES strings:
if __name__ == '__main__':
    test_smiles = [
        "O=C(OC[C@@H](O)CO)CCCCCCC/C=C(\\CCCCCCCC)/C",  # 2,3-dihydroxypropyl (Z)-10-methyloctadec-9-enoate (expected TRUE)
        "CCCCCCCC(=O)OCC(O)CO",                         # 1-monooctanoylglycerol (expected TRUE)
        "CCCCCCCCCCCCCCCCCC(=O)OC[C@H](O)CO",           # 3-stearoyl-sn-glycerol (expected TRUE)
        "O(C(=O)CCCCCCCCCCC/C=C\\C/C=C\\CCCCC)C[C@@H](O)CO",  # MG(22:2(13Z,16Z)/0:0/0:0) (expected TRUE)
        "O(C[C@@H](O)CO)C(=O)C",                        # (R)-glycerol 1-acetate (expected TRUE)
        "O(CC(O)CO)C(=O)CC",                           # Glycerol 1-propanoate (expected TRUE)
        # Some complex structures that should not be classified as monoradylglycerol:
        "O([C@@H]1[C@@H](NC(=O)C)[C@@H](O[C@@H]([C@H]1O)CO)OC[C@H]2O[C@@H](O)[C@H](NC(=O)C)[C@@H](O)[C@H]2O)",  # example sugar derivative
    ]
    for s in test_smiles:
        result, reason = is_monoradylglycerol(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n")