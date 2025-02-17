"""
Classifies: CHEBI:75659 O-acyl-L-carnitine
"""
#!/usr/bin/env python
"""
Classifies: O-acyl-L-carnitine
Definition: An O-acylcarnitine in which the carnitine component has L-configuration.
A genuine O–acyl–L–carnitine contains a carnitine fragment where a chiral carbon (with CIP “R”)
is connected (ignoring hydrogens) to exactly three heavy atoms:
  1. an oxygen that is acylated (i.e. the oxygen is in an ester bond to a carbonyl carbon),
  2. a trimethylammonium fragment (a nitrogen with formal charge +1 and three carbon neighbours),
  3. and a chain leading to a carboxylate group (e.g. –CH2–C(=O)[O-]).
If these conditions are met the molecule is classified as an O–acyl–L–carnitine.
"""

from rdkit import Chem

def is_O_acyl_L_carnitine(smiles: str):
    """
    Determines if a molecule is an O-acyl-L-carnitine based on its SMILES string.
    
    The procedure is:
      1. Parse the SMILES and assign stereochemistry.
      2. Loop through heavy atoms that are carbon and are marked as stereo centers.
      3. For each candidate:
           - verify it has a CIP label "R" (which in our depiction corresponds to L–carnitine),
           - check that (ignoring hydrogens) it is connected to exactly three heavy–atom neighbours,
           - confirm one neighbour is an oxygen that is acylated (i.e. the oxygen is bound only to
             this chiral carbon and to an acyl carbon which bears at least one double-bonded oxygen),
           - check that one neighbour is a trimethylammonium group (a positively–charged nitrogen with three carbon neighbours),
           - and that one neighbour (a carbon) is connected to an oxygen anion (as part of a carboxylate).
      4. If one candidate meets all these criteria, return True plus a descriptive reason.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is identified as an O-acyl-L-carnitine, False otherwise.
        str: Explanation for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # assign stereochemistry and compute CIP labels (important for chiral centers)
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    
    # helper: return list of heavy (atomic num >1) neighbours of an atom.
    def heavy_neighbors(atom):
        return [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
    
    # loop over all atoms that are carbons (atomic number 6) and have an assigned chiral tag
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        # require an explicit chiral tag
        if atom.GetChiralTag() == Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
            continue
        # retrieve CIP label if available
        try:
            cip = atom.GetProp('_CIPCode')
        except KeyError:
            continue
        if cip != "R":
            # We only want the candidate to have "R" CIP (i.e. L–carnitine)
            continue
        
        # Get heavy-atom neighbours (ignore hydrogens)
        nbrs = heavy_neighbors(atom)
        if len(nbrs) != 3:
            # The carnitine chiral center should be connected to exactly 3 heavy atoms
            continue

        # Initialize flags for the three groups we expect
        acyl_oxygen_ok = False
        trimethyl_ok = False
        carboxylate_ok = False

        # iterate over neighbours to check for the distinct groups
        for nbr in nbrs:
            atomic_num = nbr.GetAtomicNum()
            if atomic_num == 8:
                # Candidate for acyl oxygen.
                # Check that this oxygen is connected to exactly 2 heavy atoms:
                nbr_heavy = heavy_neighbors(nbr)
                if len(nbr_heavy) != 2:
                    continue
                # One of the neighbours must be the candidate chiral center; the other is the acyl carbon.
                other = [x for x in nbr_heavy if x.GetIdx() != atom.GetIdx()]
                if not other:
                    continue
                acyl_carb = other[0]
                if acyl_carb.GetAtomicNum() != 6:
                    continue
                # Check that the acyl carbon is part of a carbonyl: i.e.
                # at least one neighbour (other than the oxygen we came from) is an oxygen bound by a double bond.
                carbonyl_found = False
                for nn in acyl_carb.GetNeighbors():
                    if nn.GetIdx() == nbr.GetIdx():
                        continue
                    bond = mol.GetBondBetweenAtoms(acyl_carb.GetIdx(), nn.GetIdx())
                    if nn.GetAtomicNum() == 8 and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        carbonyl_found = True
                        break
                if carbonyl_found:
                    acyl_oxygen_ok = True

            elif atomic_num == 7:
                # candidate for trimethylammonium group. Check that it carries a positive charge.
                if nbr.GetFormalCharge() != 1:
                    continue
                # also, expect exactly three heavy–atom (carbon) neighbours attached to this N.
                nbrs_of_N = heavy_neighbors(nbr)
                # count only carbons among them
                methyl_count = sum(1 for x in nbrs_of_N if x.GetAtomicNum() == 6)
                if methyl_count == 3:
                    trimethyl_ok = True

            elif atomic_num == 6:
                # candidate for the chain leading to the carboxylate.
                # For our purposes we check if this carbon (likely CH2) is attached (other than to the candidate)
                # to an oxygen bearing a negative formal charge.
                for sub_nbr in nbr.GetNeighbors():
                    if sub_nbr.GetIdx() == atom.GetIdx():
                        continue
                    if sub_nbr.GetAtomicNum() == 8 and sub_nbr.GetFormalCharge() == -1:
                        carboxylate_ok = True
                        break
        if acyl_oxygen_ok and trimethyl_ok and carboxylate_ok:
            return True, ("Molecule contains a carnitine chiral center with CIP 'R', a correctly acylated "
                          "oxygen, a trimethylammonium group, and a chain leading to a carboxylate. "
                          "Thus it is classified as an O–acyl–L–carnitine.")
    
    # If we came here, then either we found a carnitine-like fragment with wrong configuration
    # or the acylation/other attachments are missing.
    return False, ("No valid O-acyl-L-carnitine motif found, or the carnitine stereocenter is not 'R' "
                   "or is missing proper acylation/trimethylammonium or carboxylate fragments.")

# Example usage (can be run as a script):
if __name__ == "__main__":
    test_smiles = [
        "CCCC(=O)O[C@H](CC([O-])=O)C[N+](C)(C)C",  # O-butanoyl-L-carnitine - expected True
        "C[N+](C)(C)C[C@@H](CC([O-])=O)OC(=O)CC(O)=O",  # O-malonyl-L-carnitine - expected True
        "O([C@@H](C[N+](C)(C)C)CC([O-])=O)C(=O)CCCCCCCCCCCC",  # CAR(13:0) false positive in previous code.
        "C[N+](C)(C)C[C@H](CC([O-])=O)OC(=O)C=C",  # O-propenoyl-D-carnitine - expected to be rejected
        "O[C@@](C([N+](CC(=O)CCC)(C)C)([2H])[2H])(CC([O-])=O)[2H]"  # Example with isotopes (butyryl-L-carnitine-d3) - expected True
    ]
    for s in test_smiles:
        result, reason = is_O_acyl_L_carnitine(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n")