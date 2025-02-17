"""
Classifies: CHEBI:75659 O-acyl-L-carnitine
"""
#!/usr/bin/env python
"""
Classifies: O-acyl-L-carnitine
Definition: An O-acylcarnitine in which the carnitine component has L-configuration.
A genuine O–acyl–L–carnitine is characterized by the presence of a carnitine fragment
where the alcohol (chiral) carbon (with CIP “R”) is directly connected to three groups:
  1. an oxygen whose only other connection is to a carbonyl carbon (i.e. it is acylated),
  2. a methylene branch that is directly attached to a trimethylammonium (N+ with three methyl groups),
  3. and a methylene branch that directly leads to a carboxylate fragment (a –COO– group).
If these conditions are met the molecule is classified as an O–acyl–L–carnitine.
"""

from rdkit import Chem

def is_O_acyl_L_carnitine(smiles: str):
    """
    Determines if a molecule is an O-acyl-L-carnitine based on its SMILES string.

    The algorithm:
      1. Parse the SMILES and assign stereochemistry.
      2. Loop over candidate chiral carbon atoms that have an assigned CIP code.
      3. For each candidate having CIP code "R", require that it has exactly three heavy-atom neighbours.
         Then, require that:
           (a) One neighbour is an oxygen that is directly esterified – i.e. it is bonded only to the chiral carbon
               and to a carbon which bears a double-bonded oxygen.
           (b) One neighbour is a carbon which (directly or via one bond) is connected to a carboxylate fragment,
               defined as a carbon bearing at least one double-bonded oxygen and bonded to an oxygen with formal charge -1.
           (c) One neighbour is a carbon that is directly connected to a nitrogen that is trimethylammonium,
               i.e. a nitrogen with formal charge +1 and three carbon neighbours.
      4. If such a chiral center is found the function returns True with an explanation.
      Otherwise returns False.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is identified as an O-acyl-L-carnitine, False otherwise.
        str: Explanation for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure stereochemical information (including CIP labels) is set.
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    
    # Helper: returns heavy (atomic number > 1) neighbors of an atom.
    def heavy_neighbors(atom):
        return [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
    
    # Helper: test whether an oxygen is “acylated” (i.e. ester oxygen) by checking that
    # its only other neighbor (aside from the chiral carbon) is a carbon that bears a double-bonded oxygen.
    def is_acylated_oxygen(o_atom, attached_to_idx):
        # o_atom should be oxygen.
        if o_atom.GetAtomicNum() != 8:
            return False
        nbrs = [n for n in heavy_neighbors(o_atom) if n.GetIdx() != attached_to_idx]
        # In an ester the oxygen should have only one other heavy neighbor.
        if len(nbrs) != 1:
            return False
        c_atom = nbrs[0]
        if c_atom.GetAtomicNum() != 6:
            return False
        # Look among carbon neighbours for a double bond (=O)
        for nbr in heavy_neighbors(c_atom):
            if nbr.GetIdx() == o_atom.GetIdx():
                continue
            bond = mol.GetBondBetweenAtoms(c_atom.GetIdx(), nbr.GetIdx())
            if nbr.GetAtomicNum() == 8 and bond is not None and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                return True
        return False

    # Helper: check for direct trimethylammonium pattern.
    def is_direct_trimethylammonium(v_atom):
        # Expect a nitrogen with formal charge +1 and at least three carbon neighbors.
        if v_atom.GetAtomicNum() == 7 and v_atom.GetFormalCharge() == 1:
            cnt = 0
            for nbr in heavy_neighbors(v_atom):
                if nbr.GetAtomicNum() == 6:
                    cnt += 1
            if cnt >= 3:
                return True
        return False

    # Helper: check a carbon for a carboxylate group.
    # A carboxylate fragment: a carbon atom attached to at least one oxygen by a double bond and also to an oxygen with -1 charge.
    def is_carboxylate_fragment(c_atom):
        if c_atom.GetAtomicNum() != 6:
            return False
        has_dcarbonyl = False
        has_neg_oxygen = False
        for nbr in heavy_neighbors(c_atom):
            bond = mol.GetBondBetweenAtoms(c_atom.GetIdx(), nbr.GetIdx())
            if nbr.GetAtomicNum() == 8:
                if bond is not None and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    has_dcarbonyl = True
                if nbr.GetAtomicNum() == 8 and nbr.GetFormalCharge() == -1:
                    has_neg_oxygen = True
        return has_dcarbonyl and has_neg_oxygen
    
    # Helper: given a candidate neighbor (usually a CH2) directly attached to the chiral center,
    # check (directly or via its own neighbor) for a carboxylate fragment.
    def has_carboxylate(candidate, parent_idx):
        # Direct: candidate itself might be the carboxylate carbon.
        if is_carboxylate_fragment(candidate):
            return True
        # Or candidate might be a methylene that leads directly to a carboxylate.
        for nbr in heavy_neighbors(candidate):
            if nbr.GetIdx() == parent_idx:
                continue
            if is_carboxylate_fragment(nbr):
                return True
        return False

    # Helper: from a candidate neighbor (expected to be carbon) check for trimethylammonium,
    # either directly or via a one-bond extension.
    def has_direct_ammonium(candidate, parent_idx):
        # Directly attached nitrogen.
        for nbr in heavy_neighbors(candidate):
            if nbr.GetIdx() == parent_idx:
                continue
            if is_direct_trimethylammonium(nbr):
                return True
        return False

    # Loop over atoms and look for a chiral carbon (atomic number 6) with an explicit CIP code.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        if atom.GetChiralTag() == Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
            continue
        try:
            cip = atom.GetProp('_CIPCode')
        except KeyError:
            continue
        # For carnitine the CIP code should be “R” (corresponding to L-carnitine).
        if cip != "R":
            continue
        
        nbrs = heavy_neighbors(atom)
        # In genuine carnitine the chiral (alcohol) carbon is bond to three heavy atoms.
        if len(nbrs) != 3:
            continue

        found_acyl = False
        found_ammonium = False
        found_carboxylate = False

        # We now expect exactly one neighbor to be the ester oxygen.
        for nbr in nbrs:
            if nbr.GetAtomicNum() == 8:
                # Check that this oxygen is directly acylated.
                if is_acylated_oxygen(nbr, atom.GetIdx()):
                    found_acyl = True
                    continue
        # For the other two neighbours, they should be carbons.
        for nbr in nbrs:
            if nbr.GetAtomicNum() == 6:
                # One branch should lead to a carboxylate group.
                if not found_carboxylate and has_carboxylate(nbr, atom.GetIdx()):
                    found_carboxylate = True
                    continue
                # The other branch should lead directly to a trimethylammonium group (one-bond away).
                if not found_ammonium and has_direct_ammonium(nbr, atom.GetIdx()):
                    found_ammonium = True
                    continue

        if found_acyl and found_ammonium and found_carboxylate:
            reason = ("Molecule contains a carnitine chiral center with CIP 'R' that is directly attached to an acylated oxygen, "
                      "a methylene leading to a carboxylate fragment, and a methylene leading directly to a trimethylammonium group. "
                      "Thus it is classified as an O–acyl–L–carnitine.")
            return True, reason

    return False, ("No valid O-acyl-L-carnitine motif found with the required direct substituents on the chiral center: "
                   "an ester oxygen (acylated), a branch leading to a carboxylate, and a branch leading to a trimethylammonium group.")

# Example usage: (When run as a script, one can test sample SMILES strings.)
if __name__ == "__main__":
    test_smiles = [
       "CCCC(=O)O[C@H](CC([O-])=O)C[N+](C)(C)C",  # O-butanoyl-L-carnitine (should be True)
       "C[N+](C)(C)C[C@@H](CC([O-])=O)OC(=O)CC(O)=O",  # O-malonyl-L-carnitine (True)
       "CCCCCCCCCCCCCCCCCC(=O)O[C@H](CC([O-])=O)C[N+](C)(C)C",  # O-octadecanoyl-L-carnitine (True)
       "CCCCCCCCCCCCCCCC(=O)O[C@H](CC([O-])=O)C[N+](C)(C)C",  # O-palmitoyl-L-carnitine (True)
       "CCC(=O)O[C@H](CC([O-])=O)C[N+](C)(C)C",  # O-propanoyl-L-carnitine (True)
       "CCCCCC\\C=C/C\\C=C/CCCCCCCC(=O)O[C@H](CC([O-])=O)C[N+](C)(C)C",  # O-linoleyl-L-carnitine (True)
       # False positive examples (should all return False):
       "O([C@@H](C[N+](C)(C)C)CC([O-])=O)C(=O)CCCCCCCCCCCC",   # CAR(13:0)
       "O([C@@H](C[N+](C)(C)C)CC([O-])=O)C(=O)CCCCCCCCCC",      # CAR(11:0)
    ]
    for s in test_smiles:
        result, reason = is_O_acyl_L_carnitine(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n")