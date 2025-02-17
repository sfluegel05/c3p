"""
Classifies: CHEBI:17088 monoacyl-sn-glycerol 3-phosphate
"""
"""
Classifies: monoacyl-sn-glycerol 3-phosphate

Definition:
  An sn-glycero-3-phosphate compound that has exactly one acyl ester
  (an acyl group replacing one –OH at either sn-1 or sn-2) while the phosphate
  group (at sn-3) remains the only additional ester substituent. Additional
  heteroatoms (for example nitrogen from choline, ethanolamine, serine, etc.) 
  are not allowed.
  
The approach is to:
  1. Parse the SMILES.
  2. Require that the molecule contains exactly one phosphorus atom.
  3. Reject molecules that contain nitrogen atoms.
  4. Count the number of acyl ester groups (using the pattern "OC(=O)[#6]") and
     require exactly one.
  5. Look for either an sn-1–acyl or sn-2–acyl glycerol–phosphate backbone.
     The SMARTS below demand three connected carbons where (for sn1) the first
     carbon bears an acyl ester and the second bears a free hydroxyl or (for sn2)
     vice‐versa. In both cases the third carbon is bound to the phosphate group.
  6. As an extra check the free –OH oxygen (atom map “oh”) is examined to ensure
     that it is not “over‐substituted” (i.e. it should have only one heavy–atom
     neighbor, its glycerol carbon).
     
If all tests are satisfied, the molecule is accepted.
"""

from rdkit import Chem

def is_monoacyl_sn_glycerol_3_phosphate(smiles: str):
    """
    Determines if a molecule is a monoacyl-sn-glycerol 3-phosphate.
    
    A monoacyl-sn-glycerol 3-phosphate (also often called lysophosphatidic acid)
    is defined as a glycerol-3-phosphate molecule in which exactly one of the
    sn-1 or sn-2 hydroxyl groups is acylated.
    
    The test requires:
      - exactly one phosphorus atom (for the phosphate)
      - no nitrogen present (to avoid e.g. choline, ethanolamine, serine, inositol headgroups)
      - exactly one acyl ester group (matching the SMARTS "OC(=O)[#6]")
      - the presence of one (and only one) glycerol-phosphate backbone with acylation 
        at sn-1 or sn-2. Two SMARTS patterns are used:
           • sn1 pattern (acyl at sn-1): "[C:1](OC(=O)[*:a])[C:2]([O:oh])[C:3]OP(=O)(O)O"
           • sn2 pattern (acyl at sn-2): "[C:1]([O:oh])[C:2](OC(=O)[*:a])[C:3]OP(=O)(O)O"
        In both patterns the free hydroxyl (atom mapped as 'oh') is explicitly labeled.
        Finally, we check that the free OH atom is indeed “free” (has only one heavy–atom neighbor).
        
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): A tuple containing a boolean decision and a reason for classification.
    """
    # 1. Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # 2. Ensure exactly one phosphorus atom is present.
    phosphorus_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if len(phosphorus_atoms) != 1:
        return False, f"Molecule must contain exactly one phosphorus atom (found {len(phosphorus_atoms)})"
    
    # 3. Reject if any nitrogen atoms are present (eliminates many unwanted headgroups).
    nitrogen_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]
    if nitrogen_atoms:
        return False, "Molecule contains nitrogen atoms suggesting an alternative headgroup"
    
    # 4. Count acyl ester bonds.
    # We define an acyl ester as an oxygen connected to a carbonyl: OC(=O)[#6]
    acyl_pattern = Chem.MolFromSmarts("OC(=O)[#6]")
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)
    if len(acyl_matches) != 1:
        return False, f"Molecule must contain exactly one acyl ester group (found {len(acyl_matches)})"
    
    # 5. Define SMARTS patterns for the glycerol-phosphate backbone.
    # Atom mapping is used so that the free hydroxyl oxygen (tagged as [O:oh])
    # can be checked later.
    # sn1: acylated at sn-1; free OH at sn-2; phosphate on sn-3.
    sn1_smarts = "[C:1](OC(=O)[*:a])[C:2]([O:oh])[C:3]OP(=O)(O)O"
    # sn2: free OH at sn-1; acylated at sn-2; phosphate on sn-3.
    sn2_smarts = "[C:1]([O:oh])[C:2](OC(=O)[*:a])[C:3]OP(=O)(O)O"
    
    patt_sn1 = Chem.MolFromSmarts(sn1_smarts)
    patt_sn2 = Chem.MolFromSmarts(sn2_smarts)
    
    matches_sn1 = mol.GetSubstructMatches(patt_sn1)
    matches_sn2 = mol.GetSubstructMatches(patt_sn2)
    total_matches = len(matches_sn1) + len(matches_sn2)
    
    # Exactly one backbone match is required.
    if total_matches == 0:
        return False, "No glycerol-3-phosphate backbone with appropriate acyl substitution found"
    if total_matches > 1:
        return False, f"Multiple ({total_matches}) backbone matches found; expected exactly one"
    
    # 6. Check that the free hydroxyl oxygen is indeed free (i.e. not further substituted).
    # Identify the match and extract the atom corresponding to [O:oh].
    backbone_match = matches_sn1[0] if matches_sn1 else matches_sn2[0]
    # In our SMARTS the free hydroxyl oxygen is not directly in the match tuple because it is not a heavy atom
    # in the main backbone pattern.
    # So we retrieve the free-O by following from the corresponding glycerol carbon.
    # For sn1, atom mapping [C:2] should have a neighbor that is an oxygen representing the free OH.
    # Similarly for sn2.
    # We attempt to locate an oxygen neighbor with atomic number 8 that is not part of an ester.
    free_oh_found = False
    for atom_idx in backbone_match:  # backbone_match contains indices for mapped carbons [C:1], [C:2], and [C:3]
        atom = mol.GetAtomWithIdx(atom_idx)
        # Look for oxygen neighbors that are not involved in a C=O (part of an acyl ester)
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 8:
                # Check that this oxygen does not have a double-bond to carbon (which would mark a carbonyl)
                is_carbonyl = False
                for sub_nbr in nbr.GetNeighbors():
                    if sub_nbr.GetAtomicNum() == 6:
                        for bond in nbr.GetBonds():
                            if bond.GetBeginAtom() == nbr or bond.GetEndAtom() == nbr:
                                if bond.GetBondTypeAsDouble() == 2:
                                    is_carbonyl = True
                                    break
                    if is_carbonyl:
                        break
                # We further demand that the oxygen itself is not bonded to any other heavy atoms
                # (free hydroxyl oxygen should have a degree of 1 when excluding implicit hydrogens).
                if not is_carbonyl and nbr.GetDegree() == 1:
                    free_oh_found = True
                    break
        if free_oh_found:
            break

    if not free_oh_found:
        return False, "Glycerol backbone free hydroxyl group is not in the expected form"
    
    return True, "Molecule has a glycerol-phosphate backbone with exactly one acyl ester group at sn-1 or sn-2"

# For testing purposes:
if __name__ == "__main__":
    # Examples of positively classified structures:
    test_smiles_list = [
        "CCCCCCCCCCCCCCCC(=O)OC[C@@H](O)COP(O)(O)=O",  # 1-icosanoyl-sn-glycero-3-phosphate
        "P(OC[C@H](O)COC(=O)CCCCCCCCCCCC)(O)(O)=O",     # 1-nonadecanoyl-glycero-3-phosphate
    ]
    for smi in test_smiles_list:
        result, reason = is_monoacyl_sn_glycerol_3_phosphate(smi)
        print(smi, result, reason)