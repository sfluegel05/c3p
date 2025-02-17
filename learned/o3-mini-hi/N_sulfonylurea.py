"""
Classifies: CHEBI:76983 N-sulfonylurea
"""
"""
Classifies: N-sulfonylurea
Definition: An N-sulfonylurea is defined as a urea moiety (–N–C(=O)–N–)
in which one (and only one) of the two urea nitrogen atoms is directly substituted 
by a sulfonyl group (–S(=O)(=O)–). This motif is common to several herbicides 
and type 2 diabetes drugs.
"""

from rdkit import Chem

def is_N_sulfonylurea(smiles: str):
    """
    Determines if a molecule is an N-sulfonylurea based on its SMILES string.
    
    Our strategy:
      1. Add hydrogens to get explicit information on the substitution.
      2. Identify candidate urea fragments by iterating over carbon atoms that
         are carbonyl carbons (i.e. have a double bond to O) and are connected to 
         exactly two nitrogen atoms.
      3. For each candidate urea fragment, examine the two urea N atoms:
         - For each, look at its substituents (other than the urea carbon) to see if 
           it is directly bonded (by a single bond) to a sulfur that itself is bound 
           to at least two oxygens by a double bond.
         - We relax requirements on hydrogen count/charge so as not to miss cases 
           where the N may be deprotonated (e.g. anionic forms).
      4. Accept the urea fragment only if exactly one of its two N atoms is sulfonylated.
    
    Args:
       smiles (str): SMILES string of the molecule.
    
    Returns:
       bool: True if at least one valid N-sulfonylurea fragment is found, False otherwise.
       str: Explanation for the outcome.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens to be able to count them properly
    mol = Chem.AddHs(mol)
    
    # Iterate over all atoms: candidate carbonyl carbons (C=O)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue  # not carbon
        
        # Look for at least one double-bond to oxygen (i.e. a C=O)
        o_double = []
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 8:
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                if bond is not None and bond.GetBondType() == Chem.BondType.DOUBLE:
                    o_double.append(nbr)
        if not o_double:
            continue  # not a carbonyl
        
        # Check connectivity: urea carbon should be bonded to exactly 2 nitrogens.
        n_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 7]
        if len(n_neighbors) != 2:
            continue
        
        # Now, candidate urea fragment found: atom is the C, and its two neighbors are the two N's.
        sulfonylated_count = 0
        sulfonylated_idx = None
        for n in n_neighbors:
            # Do not be too strict on hydrogen or formal charge – allow deprotonated N.
            valid_sulfonyl = False
            # Examine neighbors of this N excluding the urea carbon.
            for sub in n.GetNeighbors():
                if sub.GetIdx() == atom.GetIdx():
                    continue
                # Check if this neighbor is sulfur (atomic num 16) and is connected via a single bond.
                if sub.GetAtomicNum() == 16:
                    bond_NS = mol.GetBondBetweenAtoms(n.GetIdx(), sub.GetIdx())
                    if bond_NS is None or bond_NS.GetBondType() != Chem.BondType.SINGLE:
                        continue
                    # Now count how many oxygen neighbors of this sulfur are attached by a double bond.
                    dO_count = 0
                    for s_nbr in sub.GetNeighbors():
                        if s_nbr.GetAtomicNum() == 8:
                            bond_SO = mol.GetBondBetweenAtoms(sub.GetIdx(), s_nbr.GetIdx())
                            if bond_SO is not None and bond_SO.GetBondType() == Chem.BondType.DOUBLE:
                                dO_count += 1
                    if dO_count >= 2:
                        valid_sulfonyl = True
                        break  # Found a valid sulfonyl substituent attached to this N.
            if valid_sulfonyl:
                sulfonylated_count += 1
                sulfonylated_idx = n.GetIdx()
        
        # We require that exactly one of the two urea N atoms is sulfonylated.
        if sulfonylated_count == 1:
            return True, (
                f"Molecule contains a valid N-sulfonylurea motif: "
                f"urea fragment (C(=O) between N atoms) with N (atom idx {sulfonylated_idx}) "
                f"substituted by a sulfonyl group."
            )
    
    return False, "N-sulfonylurea moiety not found in the molecule."

# For testing; remove this part if using as a module.
if __name__ == "__main__":
    # Example test using glyburide
    test_smiles = "COc1ccc(Cl)cc1C(=O)NCCc1ccc(cc1)S(=O)(=O)NC(=O)NC1CCCCC1"
    result, reason = is_N_sulfonylurea(test_smiles)
    print("Result:", result)
    print("Reason:", reason)