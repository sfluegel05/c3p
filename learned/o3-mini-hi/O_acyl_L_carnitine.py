"""
Classifies: CHEBI:75659 O-acyl-L-carnitine
"""
#!/usr/bin/env python
"""
Classifies: O-acyl-L-carnitine
Definition: An O-acylcarnitine in which the carnitine component has L-configuration.
A genuine O–acyl–L–carnitine is characterized by the presence of a carnitine fragment
where a chiral carbon (with CIP “R”) is connected (ignoring hydrogens) to three groups:
  1. an oxygen that is acylated (i.e. the oxygen is in an ester bond to a carbonyl carbon),
  2. a trimethylammonium fragment (either directly or via a short alkyl linker),
  3. and a chain that leads to a carboxylate group (i.e. ultimately to a carbonyl with a negatively charged oxygen).
If these conditions are met the molecule is classified as an O–acyl–L–carnitine.
"""

from rdkit import Chem

def is_O_acyl_L_carnitine(smiles: str):
    """
    Determines if a molecule is an O-acyl-L-carnitine based on its SMILES string.
    
    The procedure is:
      1. Parse the SMILES and assign stereochemistry.
      2. Loop over heavy atoms (carbon) marked as chiral with an assigned CIP label.
      3. For each candidate that has CIP 'R' and exactly 3 heavy atom neighbours,
         check that one neighbour (directly) is an oxygen that is acylated (attached to a carbonyl carbon),
         another neighbour (directly or via a one-atom linker) contains a trimethylammonium (N+ with three C neighbours),
         and the third neighbour (directly or via a one-atom linker) leads to a carboxylate group (a carbonyl carbon with an adjacent negatively charged oxygen).
      4. If one candidate meets all these criteria, return True plus a descriptive reason.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is identified as an O-acyl-L-carnitine, False otherwise.
        str: Explanation for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    # assign stereochemistry and compute CIP labels
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    
    # helper: get heavy (atomic num >1) neighbours
    def heavy_neighbors(atom):
        return [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
    
    # helper: check for trimethylammonium pattern in an atom or in its immediate neighbour if needed
    def has_trimethylammonium(atom):
        # if atom is nitrogen, check immediately
        if atom.GetAtomicNum() == 7:
            if atom.GetFormalCharge() == 1:
                # count heavy neighbors that are carbons
                if sum(1 for x in heavy_neighbors(atom) if x.GetAtomicNum() == 6) >= 3:
                    return True
        # if atom is carbon, maybe check one bond away
        if atom.GetAtomicNum() == 6:
            for sub in heavy_neighbors(atom):
                if sub.GetAtomicNum() == 7 and sub.GetFormalCharge() == 1:
                    if sum(1 for x in heavy_neighbors(sub) if x.GetAtomicNum() == 6) >= 3:
                        return True
        return False
    
    # helper: check for chain leading to carboxylate in an atom or via a one-atom extension
    def has_carboxylate(atom):
        # if the atom is carbon, check its neighbours for a carbonyl and a negatively charged oxygen
        if atom.GetAtomicNum() == 6:
            for sub in heavy_neighbors(atom):
                if sub.GetAtomicNum() == 6:  # candidate carbonyl carbon
                    # check if sub has a double-bonded oxygen and also an oxygen with -1 charge
                    has_double_oxygen = False
                    has_neg_oxygen = False
                    for sub2 in heavy_neighbors(sub):
                        # ignore bond back to atom if present
                        if sub2.GetIdx() == atom.GetIdx():
                            continue
                        bond = mol.GetBondBetweenAtoms(sub.GetIdx(), sub2.GetIdx())
                        if sub2.GetAtomicNum() == 8 and bond is not None:
                            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                                has_double_oxygen = True
                            if sub2.GetAtomicNum() == 8 and sub2.GetFormalCharge() == -1:
                                has_neg_oxygen = True
                    if has_double_oxygen and has_neg_oxygen:
                        return True
        # also, if atom is not directly the carboxylate, try one bond extension (if atom is carbon)
        if atom.GetAtomicNum() == 6:
            for sub in heavy_neighbors(atom):
                if sub.GetAtomicNum() == 6 and sub.GetIdx() != atom.GetIdx():
                    if has_carboxylate(sub):
                        return True
        return False

    # loop over atoms that are carbons and have an explicit chiral tag 
    # (our candidate chiral center for the carnitine fragment)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        if atom.GetChiralTag() == Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
            continue
        try:
            cip = atom.GetProp('_CIPCode')
        except KeyError:
            continue
        # We require the CIP code to be "R" as in (R)-configuration corresponding to L-carnitine
        if cip != "R":
            continue
        
        nbrs = heavy_neighbors(atom)
        if len(nbrs) != 3:
            # The carnitine chiral center should be connected to exactly 3 heavy atoms.
            continue
        
        found_acyl = False
        found_ammonium = False
        found_carboxylate = False
        
        # Check each of the 3 groups:
        for nbr in nbrs:
            # 1. Acylated oxygen group should be an oxygen directly attached to the chiral center.
            if nbr.GetAtomicNum() == 8 and not found_acyl:
                nbrs2 = heavy_neighbors(nbr)
                # Expect the oxygen to be connected to the chiral center and one acyl carbon
                if len(nbrs2) == 2:
                    other_atoms = [x for x in nbrs2 if x.GetIdx() != atom.GetIdx()]
                    if other_atoms:
                        acyl_carb = other_atoms[0]
                        if acyl_carb.GetAtomicNum() == 6:
                            # check that acyl_carb is in a carbonyl (has a double-bonded oxygen)
                            has_carbonyl = False
                            for nn in heavy_neighbors(acyl_carb):
                                if nn.GetIdx() == nbr.GetIdx():
                                    continue
                                bond = mol.GetBondBetweenAtoms(acyl_carb.GetIdx(), nn.GetIdx())
                                if nn.GetAtomicNum() == 8 and bond is not None and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                                    has_carbonyl = True
                                    break
                            if has_carbonyl:
                                found_acyl = True
                                continue  # move to next neighbor

            # 2. Trimethylammonium group: either the neighbor itself or one bond extension must have this pattern.
            if not found_ammonium:
                if has_trimethylammonium(nbr):
                    found_ammonium = True
                    continue

            # 3. Carboxylate chain: either the neighbor or one bond extension should lead to a carboxylate.
            if not found_carboxylate:
                if has_carboxylate(nbr):
                    found_carboxylate = True
                    continue

        if found_acyl and found_ammonium and found_carboxylate:
            return True, ("Molecule contains a carnitine chiral center with CIP 'R' that connects to an acylated oxygen, "
                          "a (possibly remote) trimethylammonium group, and a chain leading to a carboxylate. "
                          "Thus it is classified as an O–acyl–L–carnitine.")
    
    return False, ("No valid O-acyl-L-carnitine motif found, or the carnitine chiral center in CIP 'R' "
                   "lacks the proper acylation, trimethylammonium, or carboxylate chain attachments.")

# Example usage (when run as a script, one can test with sample SMILES strings):
if __name__ == "__main__":
    test_smiles = [
        "CCCC(=O)O[C@H](CC([O-])=O)C[N+](C)(C)C",  # O-butanoyl-L-carnitine
        "C[N+](C)(C)C[C@@H](CC([O-])=O)OC(=O)CC(O)=O",  # O-malonyl-L-carnitine
        "CCCCCCCCCCCCCCCCCC(=O)O[C@H](CC([O-])=O)C[N+](C)(C)C",  # O-octadecanoyl-L-carnitine
        "CCCCCCCCCCCCCCCC(=O)O[C@H](CC([O-])=O)C[N+](C)(C)C",  # O-palmitoyl-L-carnitine
        "CCC(=O)O[C@H](CC([O-])=O)C[N+](C)(C)C",  # O-propanoyl-L-carnitine
        "CCCCCC\\C=C/C\\C=C/CCCCCCCC(=O)O[C@H](CC([O-])=O)C[N+](C)(C)C",  # O-linoleyl-L-carnitine
    ]
    for s in test_smiles:
        result, reason = is_O_acyl_L_carnitine(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n")