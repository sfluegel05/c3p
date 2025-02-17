"""
Classifies: CHEBI:50760 N-hydroxy-alpha-amino-acid
"""
"""
Classifies: N-hydroxy-alpha-amino-acid
Definition: Any amino acid (i.e. having exactly one free carboxyl group attached to an α–carbon)
in which at least one hydrogen on an amino substituent (whether on the α–nitrogen or on a side‐chain amino)
has been replaced by a hydroxy (-OH) group.

Our strategy:
  1. Parse the SMILES and add explicit H’s.
  2. Find exactly one free carboxyl group (using a SMARTS that matches both protonated and deprotonated forms).
  3. Identify the α–carbon as any carbon directly bonded to the carboxyl carbon.
  4. For every nitrogen in the molecule, check whether it carries an –OH (an oxygen neighbor that has at least one H).
  5. For each such N–OH, check whether it is connected to one of the candidate α–carbons along a short bond path (1 to 5 bonds).
     When the path length exceeds 1 bond we require that every intermediate atom is a carbon (and ideally sp³‐hybridized).
  6. If a valid N–OH is found, return True with an explanation; otherwise, return False.
  
This approach tries to exclude large molecules and side reactions that lead to false positives.
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_N_hydroxy_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-hydroxy-alpha-amino-acid based on its SMILES string.
    
    The algorithm:
      - Adds implicit hydrogens (so we can detect -OH groups properly).
      - Checks that there is exactly one carboxyl group defined as [CX3](=O)[OX2H1,OX1-].
      - Identifies candidate α–carbons as carbon atoms directly bonded to the carboxyl carbon.
      - For each nitrogen in the molecule that bears an -OH (i.e. is directly bonded to an oxygen that has a hydrogen),
        checks for a “short” path (1–5 bonds, with intermediate atoms being carbons with sp³ hybridization if more than one bond)
        connecting it to one of the candidate α–carbons.
    
    Args:
      smiles (str): SMILES string of the molecule.
    
    Returns:
      bool: True if the molecule is classified as an N-hydroxy-alpha-amino-acid, else False.
      str: Explanation text.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Add explicit hydrogens so we can check for O-H bonds reliably.
    mol = Chem.AddHs(mol)
    
    # Define a SMARTS for a free carboxyl group (protonated or deprotonated)
    carboxyl_smarts = "[CX3](=O)[OX2H1,OX1-]"
    carboxyl_pattern = Chem.MolFromSmarts(carboxyl_smarts)
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    if not carboxyl_matches:
        return False, "No free carboxyl group found; not an amino acid backbone."
    if len(carboxyl_matches) != 1:
        return False, f"Found {len(carboxyl_matches)} carboxyl groups; expected exactly one for a free amino acid."
    
    # Identify the carboxyl carbon as the first atom in the match.
    carboxyl_carbon_idx = carboxyl_matches[0][0]
    carboxyl_atom = mol.GetAtomWithIdx(carboxyl_carbon_idx)
    
    # Candidate α–carbons are carbon atoms bonded directly to the carboxyl carbon.
    candidate_alpha_idxs = []
    for nbr in carboxyl_atom.GetNeighbors():
        if nbr.GetAtomicNum() == 6:  # carbon
            candidate_alpha_idxs.append(nbr.GetIdx())
    if not candidate_alpha_idxs:
        return False, "No candidate α–carbon (carbon neighbor to carboxyl carbon) found."
    
    # Helper: Check if a given nitrogen atom has an -OH substituent.
    def has_N_hydroxy(nitrogen):
        # Loop over neighbors; we want an oxygen attached by a single bond that carries an explicit hydrogen.
        for nbr in nitrogen.GetNeighbors():
            if nbr.GetAtomicNum() == 8:  # oxygen
                bond = mol.GetBondBetweenAtoms(nitrogen.GetIdx(), nbr.GetIdx())
                if bond is None:
                    continue
                if bond.GetBondType() != rdchem.BondType.SINGLE:
                    continue
                # Check if the oxygen has any explicit hydrogen neighbor.
                for o_nbr in nbr.GetNeighbors():
                    if o_nbr.GetAtomicNum() == 1:
                        return True
        return False

    # Helper: Check that the bond path between alpha and candidate N is acceptable:
    # - path length between 1 and 5 bonds
    # - if path length >1, every intermediate atom must be a carbon
    #   (and optionally sp3-hybridized).
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
        # If directly attached, that is acceptable.
        if dist == 1:
            return True
        # For a longer path, ensure every intermediate atom is a carbon with sp3 hybridization.
        for idx in path[1:-1]:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6:
                return False
            # Check for sp3 hybridization.
            if atom.GetHybridization() != rdchem.HybridizationType.SP3:
                return False
        return True

    # Loop over all nitrogen atoms in the molecule.
    # For each, if it carries an -OH and is found (via a short acceptable bond path) from any candidate α–carbon, accept.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7:
            continue
        if not has_N_hydroxy(atom):
            continue
        n_idx = atom.GetIdx()
        # Check connection from any candidate alpha carbon to this N.
        for alpha_idx in candidate_alpha_idxs:
            if valid_n_distance(alpha_idx, n_idx):
                return True, ("Found a nitrogen with an -OH substituent (directly or via a short alkyl chain) "
                              "connected to an α–carbon adjacent to the carboxyl group.")
    
    return False, "No appropriately connected N-hydroxy nitrogen found in an amino acid backbone context."

# Example test cases (feel free to add more)
if __name__ == "__main__":
    test_smiles = [
        "O=C(O)[C@@H](NO)CCCCSC",                   # N-hydroxy-L-dihomomethionine
        "CSCCCCCCCCC(N(O)O)C(O)=O",                  # N,N-dihydroxyhexahomomethionine
        "NC(CCCN\\C(N)=N\\O)C(O)=O",                 # N(5)-[(E)-amino(hydroxyimino)methyl]ornithine
        "NC(CCCNC(N)=NO)C(O)=O",                     # N(5)-[amino(hydroxyimino)methyl]ornithine
        "CC(C)[C@H](NO)C(O)=O",                      # N-hydroxy-L-valine
        "CSCCCCCCCCC(NO)C(O)=O",                     # N-hydroxyhexahomomethionine
        "N[C@@H](CCCCNO)C(O)=O",                     # N(6)-hydroxy-L-lysine
        "O=C(O)[C@@H](N(O)O)CCCCCCCSC",              # N,N-dihydroxy-L-pentahomomethionine
        "CSCCCCCCCC(NO)C(O)=O",                      # N-hydroxypentahomomethionine
        "CSCCCCCCC(N(O)O)C(O)=O",                     # N,N-dihydroxytetrahomomethionine
        "ONCC(O)=O",                                # N-hydroxyglycine
        "ON(O)[C@@H](Cc1ccccc1)C(O)=O",              # N,N-dihydroxy-L-phenylalanine
    ]
    for smi in test_smiles:
        result, reason = is_N_hydroxy_alpha_amino_acid(smi)
        print(f"SMILES: {smi}\nClassification: {result}\nReason: {reason}\n")