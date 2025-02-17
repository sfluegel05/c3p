"""
Classifies: CHEBI:76575 monoradylglycerol
"""
"""
Classifies: monoradylglycerol
Definition: Any lipid that is glycerol bearing a single acyl, alkyl or alk-1-enyl substituent.
Approach: First, we try to locate a glycerol-like backbone (three sp3 carbons connected in a line,
each with one oxygen substituent). Then we require exactly one oxygen (the “substituted” one)
to have no hydrogen (i.e. its available valences are used to attach a carbon chain). From that oxygen,
we follow its non-backbone neighbor (which must be carbon) via a DFS that goes only over carbon atoms 
(and aborts if any ring is encountered). We then “assemble” a candidate fragment that should equal the 
entire molecule for a true monoradylglycerol.
"""

from rdkit import Chem

def is_monoradylglycerol(smiles: str):
    """
    Determines if a molecule is a monoradylglycerol based on its SMILES string.
    A monoradylglycerol is defined as a glycerol (CH2-CHOH-CH2) where one of the three -OH groups 
    is replaced by a single acyl/alkyl/alk-1-enyl chain, and the molecule contains only those atoms.
    
    Returns:
        bool: True if the entire heavy-atom set of the molecule matches the candidate fragment.
        str: A reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # A helper function to collect a contiguous chain of carbons (non-cyclic only)
    def dfs_chain(atom, blocked_ids, visited):
        # visit only sp3 or sp2 carbons; we assume fatty acyl chains are carbon-only
        if atom.GetAtomicNum() != 6:
            return set()
        if atom.IsInRing():
            # any ring encountered means we abort this branch
            return None  
        visited.add(atom.GetIdx())
        chain_set = {atom.GetIdx()}
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() in blocked_ids or nbr.GetIdx() in visited:
                continue
            if nbr.GetAtomicNum() == 6:
                result = dfs_chain(nbr, blocked_ids, visited)
                if result is None:
                    return None
                chain_set.update(result)
        return chain_set

    # Minimal chain length required (counting carbon atoms found in DFS)
    MIN_CHAIN_LENGTH = 2

    # We now iterate over potential "middle" carbons that could be the middle of a glycerol backbone.
    # In a free glycerol (or monoradylglycerol) the backbone is CH2-CHOH-CH2.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        if atom.GetHybridization().name != "SP3":
            continue
        # Candidate middle carbon should be connected (by a sigma bond) to exactly two sp3 carbons.
        nbr_carbons = [nbr for nbr in atom.GetNeighbors() 
                       if nbr.GetAtomicNum() == 6 and nbr.GetHybridization().name == "SP3"]
        if len(nbr_carbons) != 2:
            continue
        # For a linear backbone, the two neighbor carbons should not be directly bonded.
        if mol.GetBondBetweenAtoms(nbr_carbons[0].GetIdx(), nbr_carbons[1].GetIdx()):
            continue

        # Define the backbone as the three carbons: terminal1, middle, terminal2.
        backbone = [nbr_carbons[0], atom, nbr_carbons[1]]
        backbone_ids = set(a.GetIdx() for a in backbone)
        
        # For a glycerol backbone, each backbone carbon must have exactly one oxygen neighbor
        # that is not also part of the backbone. (It may have hydrogens as implicit; we check oxygen neighbors.)
        oxygen_neighbors = []  # will correspond to each backbone carbon in order
        valid_backbone = True
        for carbon in backbone:
            oxygens = [nbr for nbr in carbon.GetNeighbors() 
                       if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in backbone_ids]
            if len(oxygens) != 1:
                valid_backbone = False
                break
            oxygen_neighbors.append(oxygens[0])
        if not valid_backbone:
            continue

        # Count free -OH groups vs substituted oxygen.
        # We assume a free hydroxyl oxygen will have at least one attached hydrogen (even if implicit).
        free_OH_count = 0
        substituted_index = -1
        for idx, o_atom in enumerate(oxygen_neighbors):
            # GetTotalNumHs works well whether hydrogens are implicit or explicit.
            if o_atom.GetTotalNumHs() > 0:
                free_OH_count += 1
            else:
                substituted_index = idx
        # For monoradylglycerol, exactly two oxygens should be free and one substituted.
        if free_OH_count != 2 or substituted_index < 0:
            continue

        # Now, focus on the substituted oxygen.
        sub_o = oxygen_neighbors[substituted_index]
        # Determine which backbone carbon this oxygen is attached to.
        parent_carbon = None
        for nbr in sub_o.GetNeighbors():
            if nbr.GetIdx() in backbone_ids:
                parent_carbon = nbr
                break
        if parent_carbon is None:
            continue

        # Identify the substituent anchor: among sub_o's neighbors that are not in the backbone.
        sub_neighbors = [nbr for nbr in sub_o.GetNeighbors() if nbr.GetIdx() not in backbone_ids]
        if not sub_neighbors:
            continue
        # We require that the substituent anchor is carbon-based.
        sub_anchor = None
        for a in sub_neighbors:
            if a.GetAtomicNum() == 6:
                sub_anchor = a
                break
        if sub_anchor is None:
            continue

        # Now, traverse the substituent carbon chain.
        visited = set()
        chain_set = dfs_chain(sub_anchor, blocked_ids=backbone_ids, visited=visited)
        if chain_set is None:
            # encountered a ring in the chain -> reject candidate
            continue
        if len(chain_set) < MIN_CHAIN_LENGTH:
            continue  # chain too short to be a true acyl/alkyl chain

        # If the substituent is an acyl chain, it often has a carbonyl oxygen attached to its first carbon.
        # We add that oxygen if it is connected by a double bond.
        extra_chain_atoms = set()
        for cid in chain_set:
            ca = mol.GetAtomWithIdx(cid)
            for nbr in ca.GetNeighbors():
                # Check for oxygen connected by a double bond (indicating a carbonyl)
                bond = mol.GetBondBetweenAtoms(ca.GetIdx(), nbr.GetIdx())
                if nbr.GetAtomicNum() == 8 and bond is not None and bond.GetBondTypeAsDouble() == 2.0:
                    extra_chain_atoms.add(nbr.GetIdx())

        # Now we assemble our candidate monoradylglycerol fragment:
        candidate_atom_ids = set()
        # Add backbone carbons
        candidate_atom_ids.update(backbone_ids)
        # Add all three oxygen neighbors (free and substituted)
        for o in oxygen_neighbors:
            candidate_atom_ids.add(o.GetIdx())
        # Add all atoms from the substituent carbon chain and any extra carbonyl oxygens detected.
        candidate_atom_ids.update(chain_set)
        candidate_atom_ids.update(extra_chain_atoms)

        # Now, form the set of heavy-atom indices in the whole molecule (atomic number > 1)
        all_heavy_atom_ids = set(a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1)
        
        # To be a monoradylglycerol the candidate fragment should account for nearly all heavy atoms.
        # (If any extra heavy atoms are present beyond the candidate fragment, then the molecule is more complex.)
        if candidate_atom_ids == all_heavy_atom_ids:
            # Determine which position (1, 2 or 3) of the backbone is substituted.
            pos = substituted_index + 1
            reason = (f"Found glycerol backbone (atoms {', '.join(str(a.GetIdx()) for a in backbone)}) with substitution "
                      f"at position {pos}; substituent chain length = {len(chain_set)}, free -OH groups = 2.")
            return True, reason
        else:
            # If extra heavy atoms are present, then even though a glycerol-like fragment was found,
            # the molecule contains more than a single chain + glycerol group.
            diff = all_heavy_atom_ids - candidate_atom_ids
            reason = (f"Candidate fragment found (backbone atoms {', '.join(str(a.GetIdx()) for a in backbone)}) with substitution "
                      f"at position {substituted_index+1} but extra atoms ({len(diff)} additional heavy atoms) remain in the molecule.")
            # We continue search because a different backbone candidate might give an exact match
            # For now, we record but do not immediately return True.
            # (In our testing if any candidate gives an exact match we return True.)
            continue

    return False, "No glycerol backbone with one substituted oxygen (and two free –OH groups) forming an isolated candidate fragment was found."


# For example usage and testing:
if __name__ == '__main__':
    test_smiles = [
        "O=C(OC[C@@H](O)CO)CCCCCCC/C=C(\\CCCCCCCC)/C",  # 2,3-dihydroxypropyl (Z)-10-methyloctadec-9-enoate (expected TRUE)
        "CCCCCCCC(=O)OCC(O)CO",                         # 1-monooctanoylglycerol (expected TRUE)
        "CCCCCCCCCCCCCCCCCC(=O)OC[C@H](O)CO",           # 3-stearoyl-sn-glycerol (expected TRUE)
        "O(C(=O)CCCCCCCCCCC/C=C\\C/C=C\\CCCCC)C[C@@H](O)CO",  # MG(22:2(13Z,16Z)/0:0/0:0) (expected TRUE)
        "O(C[C@@H](O)CO)C(=O)C",                        # (R)-glycerol 1-acetate (expected TRUE)
        "O(CC(O)CO)C(=O)CC",                           # Glycerol 1-propanoate (expected TRUE)
        # A complex molecule that should not be classified as monoradylglycerol:
        "O([C@@H]1[C@@H](NC(=O)C)[C@@H](O[C@@H]([C@H]1O)CO)OC[C@H]2O[C@@H](O)[C@H](NC(=O)C)[C@@H](O)[C@H]2O)",
    ]
    for s in test_smiles:
        result, reason = is_monoradylglycerol(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n")