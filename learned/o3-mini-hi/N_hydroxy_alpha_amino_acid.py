"""
Classifies: CHEBI:50760 N-hydroxy-alpha-amino-acid
"""
"""
Classifies: N-hydroxy-alpha-amino-acid
Definition: Any amino acid (i.e. having exactly one free carboxyl group attached to an α–carbon)
in which at least one hydrogen on an amino substituent (whether on the α–nitrogen or on a side‐chain amino)
has been replaced by a hydroxy (-OH) group.
Our approach:
  1. Parse the SMILES and require exactly one carboxyl group (accepting both protonated and deprotonated forms).
  2. Identify the α–carbon as any carbon directly bonded to the carboxyl carbon.
  3. Look for any nitrogen that carries an –OH substituent.
  4. Accept if that nitrogen is “connected” to the α–carbon by a bond path containing between 1 and 5 bonds and if, 
     for paths longer than 1 bond, all intermediate atoms are carbons.
If any of these clues are missing then we reject the molecule.
"""
from rdkit import Chem

def is_N_hydroxy_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-hydroxy-alpha-amino-acid based on its SMILES string.
    
    The algorithm:
      - Checks that there is exactly one carboxyl group (C(=O)[O;H,-]).
      - Identifies the α–carbon as any carbon directly connected to the carboxyl carbon.
      - Searches for any nitrogen that carries an -OH substituent
        and is connected to an identified α–carbon by a bond path of 1 to 5 bonds
        (with the criteria that if the path length >1, every intermediate atom is a carbon).
    
    Args:
        smiles (str): SMILES representation of the molecule.
    
    Returns:
        bool: True if the molecule is classified as an N-hydroxy-alpha-amino-acid, False otherwise.
        str: Explanation of the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Look for a free carboxyl group; we accept either protonated or deprotonated forms.
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[O;H,-]")
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    if not carboxyl_matches:
        return False, "No carboxyl group found; not an amino acid backbone."
    if len(carboxyl_matches) != 1:
        return False, f"Found {len(carboxyl_matches)} carboxyl groups; expected exactly one for a free amino acid."
    
    # Identify the carboxyl carbon as the first atom in the match.
    carboxyl_carbon_idx = carboxyl_matches[0][0]
    carboxyl_atom = mol.GetAtomWithIdx(carboxyl_carbon_idx)
    
    # Identify candidate α–carbons. They must be carbon atoms directly attached to the carboxyl carbon.
    candidate_alpha_idxs = []
    for nbr in carboxyl_atom.GetNeighbors():
        if nbr.GetSymbol() == "C":
            candidate_alpha_idxs.append(nbr.GetIdx())
    if not candidate_alpha_idxs:
        return False, "No candidate α–carbon (carbon atom neighbor to carboxyl carbon) found."
    
    # Helper function: Check if a given nitrogen has an -OH substituent (i.e. a single bond to oxygen).
    def has_N_hydroxy(nitrogen):
        for o_nbr in nitrogen.GetNeighbors():
            if o_nbr.GetAtomicNum() == 8:  # oxygen
                bond = mol.GetBondBetweenAtoms(nitrogen.GetIdx(), o_nbr.GetIdx())
                if bond is not None and bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                    return True
        return False

    # Helper function: Verify whether the shortest bonded path between an α–carbon (source) and a nitrogen (target)
    # has an allowed length (1 to 5 bonds) and, if there are intermediate atoms (for paths longer than 1),
    # that all of them are carbons.
    def valid_n_distance(alpha_idx, n_idx):
        if alpha_idx == n_idx:
            return False
        try:
            path = Chem.GetShortestPath(mol, alpha_idx, n_idx)
        except Exception:
            return False
        if not path:
            return False
        dist = len(path) - 1  # number of bonds
        if dist < 1 or dist > 5:
            return False
        if dist == 1:
            return True
        # For paths longer than 1 bond, check that every intermediate atom is carbon.
        for idx in path[1:-1]:
            if mol.GetAtomWithIdx(idx).GetSymbol() != "C":
                return False
        return True

    # For each candidate α–carbon, search for any nitrogen (which is not necessarily directly attached
    # to it) that:
    #   1) Has an -OH substituent.
    #   2) Is connected to this α–carbon via a bond path meeting our distance and carbon-chain criteria.
    for alpha_idx in candidate_alpha_idxs:
        for atom in mol.GetAtoms():
            if atom.GetSymbol() == "N":
                if not has_N_hydroxy(atom):
                    continue
                n_idx = atom.GetIdx()
                if valid_n_distance(alpha_idx, n_idx):
                    return True, ("Found a nitrogen with an -OH substituent (directly or via a short carbon chain) "
                                  "that is connected to the α–carbon adjacent to the carboxyl group.")
    
    return False, "No appropriately connected N-hydroxy nitrogen found in an amino acid backbone context."

# Example testing (modify or remove as necessary)
if __name__ == "__main__":
    test_smiles = [
        "O=C(O)[C@@H](NO)CCCCSC",             # N-hydroxy-L-dihomomethionine
        "CSCCCCCCCCC(N(O)O)C(O)=O",            # N,N-dihydroxyhexahomomethionine
        "NC(CCCN\\C(N)=N\\O)C(O)=O",           # N(5)-[(E)-amino(hydroxyimino)methyl]ornithine
        "NC(CCCNC(N)=NO)C(O)=O",              # N(5)-[amino(hydroxyimino)methyl]ornithine
        "CC(C)[C@H](NO)C(O)=O",               # N-hydroxy-L-valine
        "CSCCCCCCCCC(NO)C(O)=O",              # N-hydroxyhexahomomethionine
        "N[C@@H](CCCCNO)C(O)=O",              # N(6)-hydroxy-L-lysine
        "O=C(O)[C@@H](N(O)O)CCCCCCCSC",        # N,N-dihydroxy-L-pentahomomethionine
        "CSCCCCCCCC(NO)C(O)=O",              # N-hydroxypentahomomethionine
        "CSCCCCCCC(N(O)O)C(O)=O",             # N,N-dihydroxytetrahomomethionine
    ]
    for smi in test_smiles:
        result, reason = is_N_hydroxy_alpha_amino_acid(smi)
        print(f"SMILES: {smi}\nClassification: {result}\nReason: {reason}\n")