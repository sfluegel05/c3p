"""
Classifies: CHEBI:75659 O-acyl-L-carnitine
"""
#!/usr/bin/env python
"""
Classifies: O-acyl-L-carnitine
Definition: An O-acylcarnitine in which the carnitine component has L-configuration.
A genuine O–acyl–L–carnitine is characterized by the presence of a carnitine fragment,
where a chiral carbon (with CIP “R”) is connected (ignoring hydrogens) to three groups:
  1. an oxygen that is acylated (that is, the oxygen is in an ester bond to a carbonyl carbon),
  2. a trimethylammonium fragment (either directly or via a one-atom linker),
  3. and a chain that leads (directly or via a one-atom linker) to a carboxylate group (a carbonyl bonded to an oxygen with a negative charge).
If these conditions are met the molecule is classified as an O–acyl–L–carnitine.
"""

from rdkit import Chem

def is_O_acyl_L_carnitine(smiles: str):
    """
    Determines if a molecule is an O-acyl-L-carnitine based on its SMILES string.
    
    The procedure is:
      1. Parse the SMILES and assign stereochemistry.
      2. Loop over candidate chiral carbon atoms that have an assigned CIP label.
      3. For each candidate with CIP code "R" and exactly three heavy-atom neighbors,
         check that:
           - one neighbor is an oxygen that is acylated (i.e. attached to a carbonyl carbon),
           - one neighbor (or its immediate linker) contains a trimethylammonium (N+ with three carbon neighbors),
           - and one neighbor (or its immediate linker) leads to a carboxylate fragment (a carbonyl with an adjacent negatively charged oxygen).
      4. If these criteria are met, the function returns True plus a descriptive reason.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is identified as an O-acyl-L-carnitine, False otherwise.
        str: Explanation for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # First assign stereochemistry (including CIP labels)
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    
    # Helper: returns heavy (atomic number > 1) neighbors
    def heavy_neighbors(atom):
        return [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
    
    # Helper: checks for a trimethylammonium pattern in the atom or its immediate neighbor.
    def has_trimethylammonium(atom):
        # If the atom is nitrogen, check immediately.
        if atom.GetAtomicNum() == 7:
            if atom.GetFormalCharge() == 1:
                if sum(1 for x in heavy_neighbors(atom) if x.GetAtomicNum() == 6) >= 3:
                    return True
        # If the atom is carbon, check one bond away.
        if atom.GetAtomicNum() == 6:
            for sub in heavy_neighbors(atom):
                if sub.GetAtomicNum() == 7 and sub.GetFormalCharge() == 1:
                    if sum(1 for x in heavy_neighbors(sub) if x.GetAtomicNum() == 6) >= 3:
                        return True
        return False

    # Helper: checks if an atom (expected to be carbon) is part of a carboxylate fragment.
    # A carboxylate fragment is defined as a carbon bonded to at least one double-bonded oxygen
    # and also bonded to an oxygen with a -1 formal charge.
    def is_carboxylate_fragment(atom):
        if atom.GetAtomicNum() != 6:
            return False
        has_double_oxygen = False
        has_neg_oxygen = False
        for nbr in heavy_neighbors(atom):
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
            if nbr.GetAtomicNum() == 8:
                if bond is not None and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    has_double_oxygen = True
                if nbr.GetAtomicNum() == 8 and nbr.GetFormalCharge() == -1:
                    has_neg_oxygen = True
        return has_double_oxygen and has_neg_oxygen

    # Helper: looks for a carboxylate fragment either directly on the atom or via a one-atom linker.
    def has_carboxylate(atom):
        # Direct check.
        if is_carboxylate_fragment(atom):
            return True
        # One-atom extension: check the heavy neighbors of atom.
        for nbr in heavy_neighbors(atom):
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != atom.GetIdx():
                if is_carboxylate_fragment(nbr):
                    return True
        return False

    # Loop over atoms that are carbon and have an explicit chiral tag.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        if atom.GetChiralTag() == Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
            continue
        try:
            cip = atom.GetProp('_CIPCode')
        except KeyError:
            continue
        # We require the CIP code to be "R" (corresponding to L-carnitine).
        if cip != "R":
            continue
        
        nbrs = heavy_neighbors(atom)
        # The carnitine chiral center should be connected to exactly three heavy atoms.
        if len(nbrs) != 3:
            continue
        
        found_acyl = False
        found_ammonium = False
        found_carboxylate = False
        
        # Check each of the three neighbors.
        for nbr in nbrs:
            # 1. Identify the acylated oxygen.
            if not found_acyl and nbr.GetAtomicNum() == 8:
                nbrs2 = [x for x in heavy_neighbors(nbr) if x.GetIdx() != atom.GetIdx()]
                if nbrs2 and nbrs2[0].GetAtomicNum() == 6:
                    acyl_carb = nbrs2[0]
                    # Verify that acyl_carb is in a carbonyl group (has a double-bonded oxygen).
                    for sub in heavy_neighbors(acyl_carb):
                        if sub.GetIdx() == nbr.GetIdx():
                            continue
                        bond = mol.GetBondBetweenAtoms(acyl_carb.GetIdx(), sub.GetIdx())
                        if sub.GetAtomicNum() == 8 and bond is not None and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                            found_acyl = True
                            break
            
            # 2. Check for trimethylammonium pattern in the neighbor or one bond extension.
            if not found_ammonium:
                if has_trimethylammonium(nbr):
                    found_ammonium = True
            
            # 3. Check if the neighbor or one-bond extension leads to a carboxylate.
            if not found_carboxylate:
                if has_carboxylate(nbr):
                    found_carboxylate = True
        
        if found_acyl and found_ammonium and found_carboxylate:
            return True, ("Molecule contains a carnitine chiral center with CIP 'R' that connects to an acylated oxygen, "
                          "a trimethylammonium group (directly or one bond away), and a chain leading (directly or via one-atom extension) "
                          "to a carboxylate fragment. Thus it is classified as an O–acyl–L–carnitine.")
    
    return False, ("No valid O-acyl-L-carnitine motif found or the carnitine chiral center with CIP 'R' "
                   "lacks one or more of the required acylated oxygen, trimethylammonium, or carboxylate attachments.")

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