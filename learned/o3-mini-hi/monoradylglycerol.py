"""
Classifies: CHEBI:76575 monoradylglycerol
"""
"""
Classifies: monoradylglycerol
Definition: Any lipid that is glycerol bearing a single acyl, alkyl or alk-1-enyl substituent 
at an unspecified position (i.e. exactly one of the 3 glycerol oxygen substituents is substituted and 
the other two are free –OH). Furthermore, once the glycerol fragment is “cut off” from the molecule, 
the remainder (the substituent) must be a single connected fragment with a minimum carbon chain length 
(minimum =2 for acyl; 8 for alkyl/alk-1-enyl linkages).
"""

from rdkit import Chem

def is_monoradylglycerol(smiles: str):
    """
    Determines if a molecule is a monoradylglycerol based on its SMILES string.
    
    The algorithm does the following:
      1. Adds hydrogens (so that free -OH groups are explicit).
      2. Iterates over candidate glycerol backbones: a contiguous three-carbon chain (not in a ring) 
         wherein each carbon has exactly one attached oxygen that is not part of the backbone.
      3. Among the three oxygen substituents, exactly two must be “free” (i.e. bearing at least one H)
         and one must be “substituted” (no H attached).
      4. We “remove” the backbone carbons and these three oxygen atoms from the heavy-atom set. The remainder 
         is taken as the substituent chain. We then check that:
           a. The substituted oxygen is connected to a single “anchor” atom that lies in the remainder.
           b. The remainder is a single connected fragment.
           c. The number of carbon atoms in that remainder (the substituent chain “length”) meets a minimum value 
              depending on the linkage type (if acyl, min=2, otherwise min=8).
      5. Finally, if the candidate fragment (glycerol fragment plus substituent) exactly equals the entire molecule’s
         heavy atoms, we return True together with a message that includes the position of substitution, type of linkage,
         and the number of chain carbons.
    
    If processing fails or no candidate meets the criteria then (False, <reason>) is returned.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is classified as monoradylglycerol, False otherwise.
        str: Reason for classification.
        (In case of unexpected problems, (None, None) may be returned.)
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False, "Invalid SMILES string"
        # Add explicit hydrogens to reliably count attached H on oxygens.
        mol = Chem.AddHs(mol)
    
        # Get all heavy atom indices (atomic num > 1)
        heavy_atoms = set(atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1)
    
        # Helper: Get connected component within a set of atom indices starting from seed.
        def get_fragment(seed, allowed_set):
            frag = set()
            stack = [seed]
            while stack:
                cur = stack.pop()
                if cur not in frag:
                    frag.add(cur)
                    atom = mol.GetAtomWithIdx(cur)
                    for nbr in atom.GetNeighbors():
                        nbr_idx = nbr.GetIdx()
                        if nbr_idx in allowed_set and nbr_idx not in frag:
                            stack.append(nbr_idx)
            return frag
    
        # Loop over atoms to find candidate middle carbon atoms that could form a glycerol backbone.
        for mid_atom in mol.GetAtoms():
            # Must be sp3 carbon.
            if mid_atom.GetAtomicNum() != 6 or mid_atom.GetHybridization().name != "SP3":
                continue
            # Find sp3 carbon neighbors (potential terminal glycerol carbons)
            nbr_carbons = [nbr for nbr in mid_atom.GetNeighbors() if nbr.GetAtomicNum() == 6 and nbr.GetHybridization().name == "SP3"]
            if len(nbr_carbons) != 2:
                continue
            # Ensure the two terminal carbons are NOT directly bonded (so that they form a linear chain rather than a cycle).
            if mol.GetBondBetweenAtoms(nbr_carbons[0].GetIdx(), nbr_carbons[1].GetIdx()):
                continue
            # Backbone candidate is the three atoms: terminal1, mid_atom, terminal2.
            backbone = [nbr_carbons[0], mid_atom, nbr_carbons[1]]
            backbone_ids = set(atom.GetIdx() for atom in backbone)
    
            # For each backbone carbon, require exactly one oxygen neighbor outside the backbone.
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
    
            # Count how many of the three oxygen substituents are “free” (have at least one hydrogen) versus substituted.
            free_OH_count = 0
            substituted_index = -1
            for idx, o_atom in enumerate(oxygen_neighbors):
                if o_atom.GetTotalNumHs() > 0:
                    free_OH_count += 1
                else:
                    substituted_index = idx
            # Must have exactly two free -OH groups and one substituted oxygen.
            if free_OH_count != 2 or substituted_index < 0:
                continue
    
            # For the candidate glycerol fragment, define its atom set: backbone carbons plus all three oxygen atoms.
            glycerol_atom_ids = backbone_ids.union(o.GetIdx() for o in oxygen_neighbors)
    
            # Now determine the substituent (the remainder) by taking the heavy atoms that are not in the glycerol fragment.
            remainder_ids = heavy_atoms - glycerol_atom_ids
            if not remainder_ids:
                # No substituent found.
                continue
    
            # Focus on the substituted oxygen:
            sub_o = oxygen_neighbors[substituted_index]
            # Identify the glycerol side of sub_o (should be exactly one and already in backbone).
            parent_carbon = None
            for nbr in sub_o.GetNeighbors():
                if nbr.GetIdx() in backbone_ids:
                    parent_carbon = nbr
                    break
            if parent_carbon is None:
                continue
            # Identify the substituent anchor: the neighbor of sub_o that is NOT the backbone.
            sub_neighbors = [nbr for nbr in sub_o.GetNeighbors() if nbr.GetIdx() not in backbone_ids]
            if len(sub_neighbors) != 1:
                continue
            sub_anchor = sub_neighbors[0]
            if sub_anchor.GetIdx() not in remainder_ids:
                continue
    
            # Check that the remainder (the substituent fragment) is connected.
            fragment = get_fragment(sub_anchor.GetIdx(), remainder_ids)
            if fragment != remainder_ids:
                # More than one fragment present outside glycerol.
                continue
            
            # Determine the linkage type.
            # We check if the substituent anchor (sub_anchor) has, aside from its bond to sub_o,
            # a neighbor oxygen connected by a double bond (a carbonyl) – if so, it is acyl.
            is_acyl = False
            carbonyl_o = None
            for nbr in sub_anchor.GetNeighbors():
                if nbr.GetIdx() == sub_o.GetIdx():
                    continue
                if nbr.GetAtomicNum() == 8:
                    bond = mol.GetBondBetweenAtoms(sub_anchor.GetIdx(), nbr.GetIdx())
                    if bond is not None and bond.GetBondTypeAsDouble() == 2.0:
                        is_acyl = True
                        carbonyl_o = nbr
                        break
            
            # Count number of carbon atoms in the substituent (the fragment).
            chain_carbons = sum(1 for idx in fragment if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
    
            # Establish minimum chain length depending on linkage type.
            if is_acyl:
                MIN_CHAIN_LENGTH = 2  # Allow very short acyl chains, e.g. acetate.
            else:
                MIN_CHAIN_LENGTH = 8
            if chain_carbons < MIN_CHAIN_LENGTH:
                continue
    
            # For acyl linkage, include the carbonyl oxygen if not already in the fragment.
            extra_ids = set()
            if is_acyl and carbonyl_o is not None and carbonyl_o.GetIdx() not in fragment:
                extra_ids.add(carbonyl_o.GetIdx())
    
            # Assemble the entire candidate fragment: glycerol part + substituent fragment (+ carbonyl if any).
            candidate_ids = glycerol_atom_ids.union(fragment).union(extra_ids)
    
            # Accept candidate if its heavy atoms equal the heavy atoms of the whole molecule.
            if candidate_ids == heavy_atoms:
                # Determine substitution position: based on which backbone carbon is connected to sub_o
                # (backbone order is as found in our candidate list: terminal, middle, terminal).
                positions = {backbone[0].GetIdx(): 1, backbone[1].GetIdx(): 2, backbone[2].GetIdx(): 3}
                sub_pos = positions.get(parent_carbon.GetIdx(), "unknown")
                chain_type = "acyl" if is_acyl else "alkyl/alk-1-enyl"
                reason = (f"Found glycerol backbone (atoms {sorted(list(backbone))}) with substitution at position {sub_pos} "
                          f"via {chain_type} linkage; substituent chain carbon count = {chain_carbons}; free -OH groups = 2.")
                return True, reason
        # If no candidate yielded an exact match:
        return False, "No glycerol backbone with one substituted oxygen (and two free –OH groups) forming an isolated candidate fragment was found."
    except Exception as e:
        return None, None

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
        # Some examples that previously were false positives or negatives.
        "O=C(CCC/C=C\\C\\C=C/C=C/C(C/C=C\\CCCCC)OO)OC(CO)CO",  # 12-HPETE 2-glyceryl ester (previously false negative)
    ]
    
    for s in test_smiles:
        result, reason = is_monoradylglycerol(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n")