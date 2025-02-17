"""
Classifies: CHEBI:76575 monoradylglycerol
"""
"""
Classifies: monoradylglycerol
Definition: Any lipid that is glycerol bearing a single acyl, alkyl or alk-1-enyl substituent (and nothing else).
Approach: 
  1. Find a glycerol-like backbone candidate defined as three sp3 carbons (in a linear array, with
     the two end carbons not directly bonded) each having exactly one oxygen neighbor not in the backbone.
  2. Among the three oxygen neighbors, exactly two must be “free” (i.e. have at least one implicit hydrogen)
     and one “substituted” (i.e. no H attached).
  3. The substituted oxygen should be attached to a carbon (the “anchor”). We then check the nature of the link:
     if this oxygen–carbon bond is an ester (the anchor bears a double‐bonded oxygen, i.e. a carbonyl) then we label
     the substituent as acyl (which allows shorter chains, e.g. acetate), otherwise it is assumed to be an alkyl linkage,
     and we require a minimum substituent (alkyl) chain length (here, chosen to be ≥8 carbons) in order to avoid mis‐classification.
  4. We perform a DFS from the substituent anchor (walking only over carbon atoms and aborting if any ring is encountered)
     to “collect” the chain. In addition, for an ester the carbonyl oxygen is added.
  5. Finally we “assemble” the candidate fragment (glycerol backbone carbons, its three oxygen substituents, 
     the substituent chain and, if appropriate, the carbonyl oxygen) and require that the set of heavy atoms in the candidate 
     exactly equals the heavy atoms in the molecule.
  
Due to the heuristic nature of this approach (and limitations in our choices for SMARTS/DFS) a molecule that is too “exotic”
may not be recognized. In that case (or if an unexpected error occurs) the function returns (None, None).
"""
from rdkit import Chem

def is_monoradylglycerol(smiles: str):
    """
    Determines if a molecule is a monoradylglycerol based on its SMILES string.
    
    The algorithm attempts to find a glycerol backbone candidate (three sp3 carbons in a row,
    each bearing one oxygen substituent) with exactly two free –OH groups and one substituted oxygen.
    The substituted oxygen must be attached to one carbon (the substituent anchor) from which we traverse a chain.
    If the substituent anchor has a double-bonded oxygen (ester linkage), we allow a short chain (min length 2),
    otherwise (ether linkage) we require a longer chain (min length 8). Finally, we require that the candidate fragment
    exactly equals the entire set of heavy atoms of the molecule.
    
    Args:
       smiles (str): SMILES string of the molecule.
       
    Returns:
       bool: True if the candidate fragment matches the entire molecule.
       str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Helper: DFS to traverse contiguous chain of carbons (no rings allowed)
    def dfs_chain(atom, blocked_ids, visited):
        if atom.GetAtomicNum() != 6:
            return set()
        if atom.IsInRing():
            return None  # abort branch: rings not allowed
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

    # Iterate over potential candidate middle carbons in a glycerol backbone.
    for atom in mol.GetAtoms():
        # Look for sp3 carbons only.
        if atom.GetAtomicNum() != 6 or atom.GetHybridization().name != "SP3":
            continue
        # Candidate middle carbon must have exactly two sp3 carbon neighbors.
        nbr_carbons = [nbr for nbr in atom.GetNeighbors() 
                       if nbr.GetAtomicNum() == 6 and nbr.GetHybridization().name == "SP3"]
        if len(nbr_carbons) != 2:
            continue
        # Ensure terminal carbons are not directly bonded (to force a chain, not a cycle)
        if mol.GetBondBetweenAtoms(nbr_carbons[0].GetIdx(), nbr_carbons[1].GetIdx()):
            continue

        backbone = [nbr_carbons[0], atom, nbr_carbons[1]]
        backbone_ids = set(a.GetIdx() for a in backbone)
        
        # For a glycerol backbone, each backbone carbon must have exactly one oxygen neighbor (not in backbone)
        oxygen_neighbors = []
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
        
        # Count free -OH (oxygen with at least one hydrogen) vs substituted oxygen (with no H)
        free_OH_count = 0
        substituted_index = -1
        for idx, o_atom in enumerate(oxygen_neighbors):
            # GetTotalNumHs works whether hydrogens are explicit or implicit.
            if o_atom.GetTotalNumHs() > 0:
                free_OH_count += 1
            else:
                substituted_index = idx
        # For monoradylglycerol, exactly two must be free OH and one is substituted.
        if free_OH_count != 2 or substituted_index < 0:
            continue

        # Focus on the substituted oxygen.
        sub_o = oxygen_neighbors[substituted_index]
        # Find its connection to the backbone (should be only one).
        parent_carbon = None
        for nbr in sub_o.GetNeighbors():
            if nbr.GetIdx() in backbone_ids:
                parent_carbon = nbr
                break
        if parent_carbon is None:
            continue
        
        # Identify substituent anchor: the non-backbone neighbor of the substituted oxygen.
        sub_neighbors = [nbr for nbr in sub_o.GetNeighbors() if nbr.GetIdx() not in backbone_ids]
        if len(sub_neighbors) != 1:
            continue
        sub_anchor = sub_neighbors[0]
        if sub_anchor.GetAtomicNum() != 6:
            continue

        # Determine linkage type. For an ester (acyl chain), the substituent oxygen is linked to a carbonyl.
        # In that case, sub_anchor should have at least one neighbor (other than sub_o) which is oxygen
        # and connected by a double bond.
        is_acyl = False
        for nbr in sub_anchor.GetNeighbors():
            if nbr.GetIdx() == sub_o.GetIdx():
                continue
            if nbr.GetAtomicNum() == 8:
                bond = mol.GetBondBetweenAtoms(sub_anchor.GetIdx(), nbr.GetIdx())
                if bond is not None and bond.GetBondTypeAsDouble() == 2.0:
                    is_acyl = True
                    carbonyl_o = nbr  # record the carbonyl oxygen
                    break

        # Now traverse the substituent chain starting at sub_anchor.
        visited = set()
        chain_set = dfs_chain(sub_anchor, blocked_ids=backbone_ids, visited=visited)
        if chain_set is None:
            continue  # ring encountered – reject candidate
        
        # Define minimum chain length (number of carbon atoms in substituent) depending on linkage type.
        if is_acyl:
            MIN_CHAIN_LENGTH = 2  # even acetate (2 carbons) is allowed
        else:
            MIN_CHAIN_LENGTH = 8  # reject short alkyl chains (e.g. 6-carbon chain is not acceptable here)
        if len(chain_set) < MIN_CHAIN_LENGTH:
            continue

        # For ester linkage, also include the carbonyl oxygen if not already in the chain.
        extra_atoms = set()
        if is_acyl:
            if carbonyl_o.GetIdx() not in chain_set:
                extra_atoms.add(carbonyl_o.GetIdx())

        # Assemble candidate fragment: backbone carbons, all three oxygen neighbors, substituent chain carbons and (if ester) the carbonyl O.
        candidate_atom_ids = set()
        candidate_atom_ids.update(backbone_ids)
        for o in oxygen_neighbors:
            candidate_atom_ids.add(o.GetIdx())
        candidate_atom_ids.update(chain_set)
        candidate_atom_ids.update(extra_atoms)
        
        # Get heavy-atom indices of whole molecule.
        all_heavy_atom_ids = set(a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1)
        
        # To accept candidate, the candidate fragment must account for ALL heavy atoms.
        if candidate_atom_ids == all_heavy_atom_ids:
            # Determine substitution position (1,2, or 3)
            pos = substituted_index + 1
            chain_type = "acyl" if is_acyl else "alkyl/alk-1-enyl"
            reason = (f"Found glycerol backbone (atoms {', '.join(str(a.GetIdx()) for a in backbone)}) with substitution "
                      f"at position {pos} via {chain_type} linkage; substituent chain length = {len(chain_set)}; free -OH groups = 2.")
            return True, reason
        # If candidate fragment does not equal full molecule, try next candidate.
    return False, "No glycerol backbone with one substituted oxygen (and two free –OH groups) forming an isolated candidate fragment was found."

# Example usage:
if __name__ == '__main__':
    test_smiles = [
        # True positives (should be classified as monoradylglycerol)
        "O=C(OC[C@@H](O)CO)CCCCCCC/C=C(\\CCCCCCCC)/C",  # 2,3-dihydroxypropyl (Z)-10-methyloctadec-9-enoate
        "CCCCCCCC(=O)OCC(O)CO",                         # 1-monooctanoylglycerol
        "CCCCCCCCCCCCCCCCCC(=O)OC[C@H](O)CO",           # 3-stearoyl-sn-glycerol
        "O(C(=O)CCCCCCCCCCC/C=C\\C/C=C\\CCCCC)C[C@@H](O)CO",  # MG(22:2(13Z,16Z)/0:0/0:0)
        "O(C[C@@H](O)CO)C(=O)C",                        # (R)-glycerol 1-acetate
        "CCCCCCCCCCCCCCCCCOC[C@@H](O)CO",               # 1-O-octadecyl-sn-glycerol
        # False positives (should not be classified as monoradylglycerol)
        "O(CCCCCC)C(CO)CO",                             # 2-Hexylglycerol (linkage is ether but chain too short: 6 carbons)
        # False negatives (should be monoradylglycerol but might be missed if backbone not recognized)
        "O=C(CCC/C=C\\C\\C=C/C=C/C(C/C=C\\CCCCC)OO)OC(CO)CO",  # 12-HPETE 2-glyceryl ester
    ]
    
    for s in test_smiles:
        result, reason = is_monoradylglycerol(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n")