"""
Classifies: CHEBI:50760 N-hydroxy-alpha-amino-acid
"""
"""
Classifies: N-hydroxy-alpha-amino-acid
Definition: Any amino acid (i.e. having exactly one free carboxyl group attached to an α–carbon)
in which at least one hydrogen on an amino substituent (whether on the α–nitrogen or on a side‐chain amino)
has been replaced by a hydroxy (-OH) group.

This implementation attempts to improve on the previous approach by:
  1. Rejecting molecules with a molecular weight above 500 Da (most free amino acids are smaller).
  2. Identifying the free carboxyl group and its attached α–carbon.
  3. Searching for any nitrogen that bears an -OH (i.e. an oxygen directly bound that has at least one hydrogen)
     and then checking if that nitrogen is reachable from one of the candidate α–carbons via a short bond path (1–5 bonds).
     (In contrast to before, we no longer require that all intermediate atoms be carbons with sp³-hybridization.)
"""

from rdkit import Chem
from rdkit.Chem import rdchem, rdMolDescriptors

def is_N_hydroxy_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is classified as an N-hydroxy-alpha-amino-acid.

    Algorithm:
      1. Parse the SMILES and add H’s.
      2. Reject molecules with molecular weight > 500 Da.
      3. Find exactly one free carboxyl group defined by SMARTS "[CX3](=O)[OX2H1,OX1-]".
      4. Identify candidate α–carbons as carbon atoms directly bonded to the carboxyl carbon.
      5. For every nitrogen that has an -OH substituent (i.e. an oxygen neighbor with at least one hydrogen attached)
         check if there is a short path (1–5 bonds) from one of the candidate α–carbons.
         (Note: we relax the requirement about the composition of intermediate atoms to allow e.g. a nitrogen in a side chain.)
      6. If such a connection is found, return True with an explanation; else return False.

    Args:
      smiles (str): SMILES string of the molecule.
    
    Returns:
      bool: True if the molecule is classified as an N-hydroxy-alpha-amino-acid, else False.
      str: Explanation text.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Add explicit hydrogens so -OH groups are represented correctly.
    mol = Chem.AddHs(mol)
    
    # Check molecular weight: many amino acids have MW < 500 Da.
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if mw > 500:
        return False, f"Molecular weight ({mw:.1f} Da) too high to be a simple amino acid."
    
    # Define a pattern for a free carboxyl group: the carbon and its oxygen(s).
    carboxyl_smarts = "[CX3](=O)[OX2H1,OX1-]"
    carboxyl_pattern = Chem.MolFromSmarts(carboxyl_smarts)
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    if not carboxyl_matches:
        return False, "No free carboxyl group found; not an amino acid backbone."
    if len(carboxyl_matches) != 1:
        return False, f"Found {len(carboxyl_matches)} free carboxyl groups; expected exactly one for a free amino acid."
    
    # Determine the carboxyl carbon from the match – use the first atom in the SMARTS match.
    carboxyl_carbon_idx = carboxyl_matches[0][0]
    carboxyl_atom = mol.GetAtomWithIdx(carboxyl_carbon_idx)
    
    # Candidate α–carbons: carbon(s) directly bonded to the carboxyl carbon.
    candidate_alpha_idxs = []
    for nbr in carboxyl_atom.GetNeighbors():
        if nbr.GetAtomicNum() == 6:  # carbon
            candidate_alpha_idxs.append(nbr.GetIdx())
    if not candidate_alpha_idxs:
        return False, "No candidate α–carbon (carbon neighbor to carboxyl carbon) found."
    
    # Helper: Check if a given nitrogen atom has an -OH substituent.
    def has_N_hydroxy(nitrogen):
        # Look for an oxygen bound by a single bond that carries an explicit hydrogen.
        for nbr in nitrogen.GetNeighbors():
            if nbr.GetAtomicNum() == 8:  # oxygen
                bond = mol.GetBondBetweenAtoms(nitrogen.GetIdx(), nbr.GetIdx())
                if not bond or bond.GetBondType() != rdchem.BondType.SINGLE:
                    continue
                # Check for at least one hydrogen attached to this oxygen.
                for subnbr in nbr.GetNeighbors():
                    if subnbr.GetAtomicNum() == 1:
                        return True
        return False

    # Helper: Check if there exists a short path (1-5 bonds) between an alpha carbon and a given nitrogen.
    def valid_path(alpha_idx, n_idx):
        # Get the shortest path (list of atom indices).
        try:
            path = Chem.GetShortestPath(mol, alpha_idx, n_idx)
        except Exception:
            return False
        if not path:
            return False
        path_length = len(path) - 1  # number of bonds
        return 1 <= path_length <= 5

    # Loop over nitrogen atoms and look for ones with an -OH substituent.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7:
            continue
        if not has_N_hydroxy(atom):
            continue
        n_idx = atom.GetIdx()
        # Check if this nitrogen is "close" (short bond path) to any candidate alpha carbon.
        for alpha_idx in candidate_alpha_idxs:
            if valid_path(alpha_idx, n_idx):
                return True, ("Found a nitrogen with an -OH substituent (directly or via a short alkyl chain) "
                              "connected (within 1-5 bonds) to an α–carbon adjacent to the carboxyl group.")
    
    # If no appropriate N-hydroxy is found.
    return False, "No appropriately connected N-hydroxy nitrogen found in an amino acid backbone context."

# Example test cases (feel free to add more)
if __name__ == "__main__":
    test_smiles = [
        "O=C(O)[C@@H](NO)CCCCSC",                   # N-hydroxy-L-dihomomethionine (TP)
        "CSCCCCCCCCC(N(O)O)C(O)=O",                  # N,N-dihydroxyhexahomomethionine (TP)
        "NC(CCCN\\C(N)=N\\O)C(O)=O",                 # N(5)-[(E)-amino(hydroxyimino)methyl]ornithine (FN before)
        "NC(CCCNC(N)=NO)C(O)=O",                     # N(5)-[amino(hydroxyimino)methyl]ornithine (FN before)
        "CC(C)[C@H](NO)C(O)=O",                      # N-hydroxy-L-valine (TP)
        "CSCCCCCCCCC(NO)C(O)=O",                     # N-hydroxyhexahomomethionine (TP)
        "N[C@@H](CCCCNO)C(O)=O",                     # N(6)-hydroxy-L-lysine (TP)
        "O=C(O)[C@@H](N(O)O)CCCCCCCSC",              # N,N-dihydroxy-L-pentahomomethionine (TP)
        "CSCCCCCCCC(NO)C(O)=O",                      # N-hydroxypentahomomethionine (TP)
        "CSCCCCCCC(N(O)O)C(O)=O",                     # N,N-dihydroxytetrahomomethionine (TP)
        "ONCC(O)=O",                                # N-hydroxyglycine (TP)
        "ON(O)[C@@H](Cc1ccccc1)C(O)=O",              # N,N-dihydroxy-L-phenylalanine (TP)
    ]
    for smi in test_smiles:
        result, reason = is_N_hydroxy_alpha_amino_acid(smi)
        print(f"SMILES: {smi}\nClassification: {result}\nReason: {reason}\n")