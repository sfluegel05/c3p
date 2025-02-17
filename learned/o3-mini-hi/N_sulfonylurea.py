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
      1. Add hydrogens so that substitution details are explicit.
      2. Use a SMARTS query to identify urea fragments of the form: [NX3]C(=O)[NX3]. 
         (This finds candidate –N–C(=O)–N– fragments.)
      3. For each candidate fragment, examine the two urea nitrogen atoms.
         For each of these N atoms, check (ignoring the bond back to the urea carbon)
         whether they are substituted by a sulfonyl group. In our approach a sulfonyl 
         substituent is defined as a sulfur atom (atomic number 16) attached via a single 
         bond that carries at least two oxygen atoms via double bonds.
      4. Accept the candidate if exactly one of the two urea nitrogens is sulfonylated. 
         (This implements the “one (and only one) of the two urea Ns” rule.)
    
    Args:
       smiles (str): SMILES string of the molecule.
       
    Returns:
       bool: True if at least one valid N-sulfonylurea fragment is found, False otherwise.
       str: Explanation for the outcome.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens to help in counting substituents
    mol = Chem.AddHs(mol)
    
    # Define a SMARTS pattern to detect the urea fragment: N-C(=O)-N.
    urea_smarts = "[NX3;!$([NX3]-C(=O)-[NX3])]"  # base pattern (we will extract the three atoms)
    # Here we use an explicit pattern for the fragment.
    # Instead, we will create the urea fragment like this:
    urea_pattern = Chem.MolFromSmarts("[NX3]C(=O)[NX3]")
    if urea_pattern is None:
        return False, "Error in urea SMARTS pattern"
    
    matches = mol.GetSubstructMatches(urea_pattern)
    if not matches:
        return False, "No urea fragment (N-C(=O)-N) found"
    
    # Helper routine: check if a given nitrogen atom (provided its idx and the index of the urea carbon)
    # is sulfonylated; i.e. has a substituent (apart from the urea carbon) which is a sulfur atom bound 
    # via a single bond that itself has at least two double-bonded oxygens.
    def is_n_sulfonylated(n_atom, urea_carbon_idx):
        for nb in n_atom.GetNeighbors():
            if nb.GetIdx() == urea_carbon_idx:
                continue  # skip the urea carbon
            # Look for a sulfur neighbor
            if nb.GetAtomicNum() == 16:
                bond = mol.GetBondBetweenAtoms(n_atom.GetIdx(), nb.GetIdx())
                if bond is None or bond.GetBondType() != Chem.BondType.SINGLE:
                    continue
                # Count number of oxygens attached to this sulfur with double bonds.
                dO_count = 0
                for s_nbr in nb.GetNeighbors():
                    # Do not count the nitrogen already attached.
                    if s_nbr.GetIdx() == n_atom.GetIdx():
                        continue
                    if s_nbr.GetAtomicNum() == 8:
                        bond_SO = mol.GetBondBetweenAtoms(nb.GetIdx(), s_nbr.GetIdx())
                        if bond_SO is not None and bond_SO.GetBondType() == Chem.BondType.DOUBLE:
                            dO_count += 1
                if dO_count >= 2:
                    # Found a valid sulfonyl substituent on this N.
                    return True
        return False

    # Now iterate over each urea fragment match.
    # Each match is a tuple: (N1, C, N2) where C is the carbonyl carbon.
    for match in matches:
        if len(match) != 3:
            continue  # should be exactly 3 atoms in the match
        n1_idx, c_idx, n2_idx = match
        n1 = mol.GetAtomWithIdx(n1_idx)
        n2 = mol.GetAtomWithIdx(n2_idx)
        
        # Verify that the central atom is a carbon with at least one double-bond to O.
        carbon = mol.GetAtomWithIdx(c_idx)
        valid_carbonyl = False
        for nb in carbon.GetNeighbors():
            if nb.GetAtomicNum() == 8:
                bond = mol.GetBondBetweenAtoms(c_idx, nb.GetIdx())
                if bond is not None and bond.GetBondType() == Chem.BondType.DOUBLE:
                    valid_carbonyl = True
                    break
        if not valid_carbonyl:
            continue  # not a proper carbonyl carbon
        
        sulfonylated_info = []  # list of (N_idx, bool)
        # For each urea nitrogen, check sulfonylation.
        for n_atom in (n1, n2):
            sulfonylated = is_n_sulfonylated(n_atom, c_idx)
            sulfonylated_info.append( (n_atom.GetIdx(), sulfonylated) )
        
        # Count how many nitrogens are sulfonylated.
        count_sulfonyl = sum(1 for idx, flag in sulfonylated_info if flag)
        if count_sulfonyl == 1:
            # Found a valid N-sulfonylurea fragment.
            sulfonylated_atom_idx = [idx for idx, flag in sulfonylated_info if flag][0]
            return True, (
                f"Molecule contains a valid N-sulfonylurea motif: "
                f"urea fragment (C(=O) between N atoms) with N (atom idx {sulfonylated_atom_idx}) "
                f"substituted by a sulfonyl group."
            )
    
    return False, "No valid N-sulfonylurea motif found in the molecule."

# For testing; remove the section below if using this as an importable module.
if __name__ == "__main__":
    # Example test using glyburide (a known N-sulfonylurea)
    test_smiles = "COc1ccc(Cl)cc1C(=O)NCCc1ccc(cc1)S(=O)(=O)NC(=O)NC1CCCCC1"
    result, reason = is_N_sulfonylurea(test_smiles)
    print("Result:", result)
    print("Reason:", reason)