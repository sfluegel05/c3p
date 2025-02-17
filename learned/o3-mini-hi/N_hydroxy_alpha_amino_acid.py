"""
Classifies: CHEBI:50760 N-hydroxy-alpha-amino-acid
"""
"""
Classifies: N-hydroxy-alpha-amino-acid
Definition: Any amino acid (i.e. a molecule having a single carboxyl group attached to 
an α–carbon) in which at least one of the hydrogens on any amino substituent (whether on the α–nitrogen or on a side-chain amino)
has been replaced by a hydroxy (-OH) group.
Our strategy is to:
  1. Locate the (free) carboxyl group using a SMARTS that matches both protonated and deprotonated forms.
  2. Identify the α–carbon as a non-oxygen atom directly attached to the carboxyl carbon.
  3. For that candidate α–carbon, look for a nitrogen (either directly attached, or along a short chain – typical of Lys/Orn side chains)
     that carries an -OH substituent. We accept a bond distance of 1 (direct attachment) or 4–5 bonds (side-chain amino groups) 
     provided that the atoms in-between are all carbons.
  4. To avoid many false positives from peptide or complex systems we require that exactly one carboxyl group is present.
  
If any of these clues are missing then we conclude the molecule is not a N-hydroxy-α–amino‐acid.
"""
from rdkit import Chem

def is_N_hydroxy_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-hydroxy-alpha-amino-acid based on its SMILES string.
    
    Our approach:
      1. Parse the SMILES string.
      2. Look for exactly one carboxyl group (accepting either –OH or –O– terminal).
      3. Identify the candidate α–carbon as a non-oxygen neighbor of that carboxyl carbon.
      4. For each candidate α–carbon, search for a nitrogen that is substituted with at least one –OH.
         We check first if such a nitrogen is directly attached (distance=1).
         If not, then we check if there is such a nitrogen at a distance of 4 or 5 bonds (typical for Lys/Orn side‐chain amino groups)
         ensuring that all intermediate atoms (if any) along the shortest path are carbons.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as an N-hydroxy-alpha-amino-acid, False otherwise.
        str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."

    # Define carboxyl SMARTS accepting both -OH and -O- forms.
    # This pattern will match a carbon double-bonded to an oxygen then bonded to an oxygen that is either neutral or negatively charged.
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[O;H,-]")
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    
    if not carboxyl_matches:
        return False, "No carboxyl group found; not a standard amino acid backbone."
    if len(carboxyl_matches) != 1:
        return False, f"Found {len(carboxyl_matches)} carboxyl groups; expected exactly one for a free amino acid."

    # Use the first match: note that in the pattern "C(=O)[O;H,-]", the first atom (index 0) is the carboxyl carbon.
    carboxyl_carbon_idx = carboxyl_matches[0][0]
    carboxyl_atom = mol.GetAtomWithIdx(carboxyl_carbon_idx)
    
    # Identify candidate α–carbon: a neighbor of the carboxyl carbon that is not oxygen.
    candidate_alpha_idxs = []
    for nbr in carboxyl_atom.GetNeighbors():
        if nbr.GetSymbol() != "O":
            candidate_alpha_idxs.append(nbr.GetIdx())
    if not candidate_alpha_idxs:
        return False, "No candidate α–carbon found (no non-oxygen neighbor to carboxyl carbon)."
    
    # Helper function to check if a given nitrogen has an -OH substituent (via a single bond).
    def has_N_hydroxy(nitrogen):
        for o_nbr in nitrogen.GetNeighbors():
            # Check for oxygen attached by a SINGLE bond.
            if o_nbr.GetAtomicNum() == 8:
                bond = mol.GetBondBetweenAtoms(nitrogen.GetIdx(), o_nbr.GetIdx())
                if bond and bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                    return True
        return False

    # Helper function: Given an alpha carbon and a target nitrogen, check that the shortest path from the α–carbon
    # to the nitrogen has the expected length (1, 4, or 5 bonds) and if longer than 1, then all intermediate atoms are carbons.
    def valid_n_distance(alpha_idx, n_idx):
        path = Chem.GetShortestPath(mol, alpha_idx, n_idx)
        if not path:
            return False
        # The number of bonds is len(path)-1.
        dist = len(path) - 1
        # We directly accept a nitrogen directly attached to the α–carbon.
        if dist == 1:
            return True
        # For side-chain amino groups in e.g. ornithine (distance = 4) or lysine (distance = 5):
        if dist in [4, 5]:
            # Check that all atoms in between (excluding endpoints) are carbons.
            for idx in path[1:-1]:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetSymbol() != "C":
                    return False
            return True
        return False

    # For each candidate α–carbon, search among the nitrogen atoms in the molecule.
    for alpha_idx in candidate_alpha_idxs:
        alpha_atom = mol.GetAtomWithIdx(alpha_idx)
        # Consider every atom that is nitrogen that carries an -OH substituent.
        for atom in mol.GetAtoms():
            if atom.GetSymbol() == "N" and has_N_hydroxy(atom):
                # Now, check the connectivity: the nitrogen should be “backbone-associated”
                # either being directly attached (distance 1) to the α–carbon or at a distance typical for a side-chain amino group.
                if valid_n_distance(alpha_idx, atom.GetIdx()):
                    return True, "Found amino nitrogen (directly or via a short side chain) with an –OH substituent, attached to an α–carbon with a carboxyl group."
    return False, "No appropriately connected N-hydroxy group found in an amino acid backbone context."

# Example testing (this block is optional and may be removed if only the function is desired):
if __name__ == "__main__":
    test_smiles = [
        # True positives
        "O=C(O)[C@@H](NO)CCCCSC",         # N-hydroxy-L-dihomomethionine (α–nitrogen directly N-hydroxy)
        "CSCCCCCCCCC(N(O)O)C(O)=O",        # N,N-dihydroxyhexahomomethionine (α–nitrogen directly N-hydroxy)
        "CC(C)[C@H](NO)C(O)=O",            # N-hydroxy-L-valine (α–nitrogen directly N-hydroxy)
        "CSCCCCCCCCC(NO)C(O)=O",           # N-hydroxyhexahomomethionine (α–nitrogen directly N-hydroxy)
        "O=C(O)[C@@H](N(O)O)CCCCCCCSC",    # N,N-dihydroxy-L-pentahomomethionine (α–nitrogen directly N-hydroxy)
        "CSCCCCCCCC(NO)C(O)=O",            # N-hydroxypentahomomethionine (α–nitrogen directly N-hydroxy)
        "CSCCCCCCC(N(O)O)C(O)=O",          # N,N-dihydroxytetrahomomethionine (α–nitrogen directly N-hydroxy)
        "O=C(O)[C@@H](N(O)O)CCCCSC",       # N,N-dihydroxy-L-dihomomethionine (α–nitrogen directly N-hydroxy)
        "CSCCCCCCCC(N(O)O)C(O)=O",         # N,N-dihydroxypentahomomethionine (α–nitrogen directly N-hydroxy)
        "O=C(O)[C@@H](N(O)O)CCCCCCCCSC",   # N,N-dihydroxy-L-hexahomomethionine (α–nitrogen directly N-hydroxy)
        "N1([C@@H](CCCC1)C(=O)O)O",        # N-hydroxy-L-pipecolic acid (cyclic but with direct N–OH on α–nitrogen)
        "CSCCCCCC(NO)C(O)=O",              # N-hydroxytrihomomethionine (α–nitrogen directly N-hydroxy)
        "O=C(O)[C@@H](NO)CCCCCCSC",        # N-hydroxy-L-tetrahomomethionine (α–nitrogen directly N-hydroxy)
        "ONCC(O)=O",                      # N-hydroxyglycine (α–nitrogen directly N-hydroxy)
        "ON(O)[C@@H](Cc1ccccc1)C(O)=O",    # N,N-dihydroxy-L-phenylalanine (α–nitrogen directly N-hydroxy)
        "O=C(O)[C@@H](NO)CCCCCSC",         # N-hydroxy-L-trihomomethionine (α–nitrogen directly N-hydroxy)
        "CC(C)[C@H](N(O)O)C(O)=O",         # N,N-dihydroxy-L-valine (α–nitrogen directly N-hydroxy)
        "CSCCCCCCC(NO)C(O)=O",             # N-hydroxytetrahomomethionine (α–nitrogen directly N-hydroxy)
        "CSCCCCCC(N(O)O)C(O)=O",           # N,N-dihydroxytrihomomethionine (α–nitrogen directly N-hydroxy)
        "CSCCCCC(N(O)O)C(O)=O",            # N,N-dihydroxydihomomethionine (α–nitrogen directly N-hydroxy)
        "CC[C@H](C)[C@H](N(O)O)C(O)=O",     # N,N-dihydroxy-L-isoleucine (α–nitrogen directly N-hydroxy)
        "ON(O)[C@@H](Cc1ccc(O)cc1)C(O)=O",  # N,N-dihydroxy-L-tyrosine (α–nitrogen directly N-hydroxy)
        "ON(O)[C@@H](Cc1c[nH]c2ccccc12)C(O)=O",  # N,N-dihydroxy-L-tryptophan (α–nitrogen directly N-hydroxy)
        # Examples that were missed previously (false negatives now accepted):
        "N[C@@H](CCCCNO)C(O)=O",           # N(6)-hydroxy-L-lysine; the N–OH is on a side-chain amino group.
        "NC(CCCN\\C(N)=N\\O)C(O)=O",       # N(5)-[(E)-amino(hydroxyimino)methyl]ornithine (side-chain amino N–OH at distance 4)
        "NC(CCCNC(N)=NO)C(O)=O",           # N(5)-[amino(hydroxyimino)methyl]ornithine (side-chain amino N–OH at distance 4)
        "N[C@@H](CCCNC(=N)NO)C(O)=O",       # N(5)-[(hydroxyamino)(imino)methyl]ornithine (side-chain amino N–OH at distance 4)
        "N[C@@H](CCCNC(=N)NO)C(O)=O",       # N(5)-[(hydroxyamino)(imino)methyl]-L-ornithine (side-chain amino N–OH at distance 4)
        "C(\\N)(=N/O)/NCCC[C@H](N)C(=O)O"   # N(5)-[(Z)-amino(hydroxyimino)methyl]-L-ornithine (side-chain amino N–OH)
    ]
    
    for smi in test_smiles:
        result, reason = is_N_hydroxy_alpha_amino_acid(smi)
        print(f"SMILES: {smi}\nClassification: {result}\nReason: {reason}\n")