"""
Classifies: CHEBI:50760 N-hydroxy-alpha-amino-acid
"""
"""
Classifies: N-hydroxy-alpha-amino-acid
Definition: Any amino acid (i.e. having a single carboxyl group attached to an α–carbon)
in which at least one hydrogen on an amino substituent (whether on the α–nitrogen or on a side-chain amino)
has been replaced by a hydroxy (-OH) group.
Our approach:
  1. Parse the SMILES string.
  2. Require exactly one carboxyl group (accept both protonated and deprotonated forms).
  3. Identify the α–carbon as a carbon (C) directly attached to the carboxyl carbon.
  4. Search for a nitrogen atom bearing an -OH substituent either directly attached to that α–carbon (distance=1)
     or attached via a chain of carbons (distance 4 or 5 bonds); in the latter case ensure intermediate atoms are all carbons.
  
If any of these clues are missing then the molecule is not classified as a N-hydroxy-alpha-amino-acid.
"""
from rdkit import Chem

def is_N_hydroxy_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-hydroxy-alpha-amino-acid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a N-hydroxy-alpha-amino-acid, False otherwise.
        str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Define a SMARTS for a carboxyl group that accepts both protonated and deprotonated forms.
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[O;H,-]")
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    if not carboxyl_matches:
        return False, "No carboxyl group found; not an amino acid backbone."
    if len(carboxyl_matches) != 1:
        return False, f"Found {len(carboxyl_matches)} carboxyl groups; expected exactly one for a free amino acid."
    
    # In our SMARTS, the first atom (index 0) is the carboxyl carbon.
    carboxyl_carbon_idx = carboxyl_matches[0][0]
    carboxyl_atom = mol.GetAtomWithIdx(carboxyl_carbon_idx)
    
    # Identify candidate α–carbon: must be a carbon (C) directly attached to the carboxyl carbon.
    candidate_alpha_idxs = []
    for nbr in carboxyl_atom.GetNeighbors():
        if nbr.GetSymbol() == "C":
            candidate_alpha_idxs.append(nbr.GetIdx())
    if not candidate_alpha_idxs:
        return False, "No candidate α–carbon (carbon atom neighbor to carboxyl carbon) found."
    
    # Helper function: Check if a given nitrogen has an -OH substituent (attached by a single bond).
    def has_N_hydroxy(nitrogen):
        for o_nbr in nitrogen.GetNeighbors():
            if o_nbr.GetAtomicNum() == 8:  # oxygen
                bond = mol.GetBondBetweenAtoms(nitrogen.GetIdx(), o_nbr.GetIdx())
                if bond and bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                    # Optionally, one might check that the oxygen is not double-bonded to any carbon
                    return True
        return False

    # Helper function: Given alpha atom (source) and a nitrogen (target), confirm that the shortest
    # path length is acceptable. Accept distance 1 (direct bond) or 4/5 (side-chain) provided that
    # all intermediate atoms (if any) are carbons.
    def valid_n_distance(alpha_idx, n_idx):
        # Avoid the degenerate case where source and target are the same.
        if alpha_idx == n_idx:
            return False
        try:
            path = Chem.GetShortestPath(mol, alpha_idx, n_idx)
        except Exception as e:
            return False
        if not path:
            return False
        # The number of bonds is one less than the number of atoms in the path.
        dist = len(path) - 1
        if dist == 1:
            return True
        if dist in [4, 5]:
            # Check that all intermediate atoms (excluding endpoints) are carbons.
            for idx in path[1:-1]:
                if mol.GetAtomWithIdx(idx).GetSymbol() != "C":
                    return False
            return True
        return False

    # For each candidate α–carbon, search for a nitrogen (other than the α–carbon itself) that:
    #   1. Has at least one -OH substituent.
    #   2. Is directly attached to the α–carbon OR attached via a chain of carbons (4 or 5 bonds away).
    for alpha_idx in candidate_alpha_idxs:
        for atom in mol.GetAtoms():
            if atom.GetSymbol() == "N":
                n_idx = atom.GetIdx()
                if not has_N_hydroxy(atom):
                    continue
                # Check connectivity only if the nitrogen is not the candidate α–carbon.
                if valid_n_distance(alpha_idx, n_idx):
                    return True, ("Found a nitrogen with an -OH substituent (either directly attached "
                                  "or via a short carbon chain) connected to the α–carbon adjacent to the carboxyl group.")
    
    return False, "No appropriately connected N-hydroxy nitrogen found in an amino acid backbone context."


# Example testing (remove or modify as needed)
if __name__ == "__main__":
    test_smiles = [
        "O=C(O)[C@@H](NO)CCCCSC",           # N-hydroxy-L-dihomomethionine
        "CSCCCCCCCCC(N(O)O)C(O)=O",          # N,N-dihydroxyhexahomomethionine
        "NC(CCCN\\C(N)=N\\O)C(O)=O",         # N(5)-[(E)-amino(hydroxyimino)methyl]ornithine
        "OC(=O)NO",                        # Problematic example that previously caused an invariant violation
        "ONCC(O)=O",                       # N-hydroxyglycine
        "N[C@@H](CCCNC(=N)NO)C(O)=O",        # N(5)-[(hydroxyamino)(imino)methyl]-L-ornithine
        "CC(C)[C@H](N(O)O)C(O)=O",          # N,N-dihydroxy-L-valine
    ]
    for smi in test_smiles:
        result, reason = is_N_hydroxy_alpha_amino_acid(smi)
        print(f"SMILES: {smi}\nClassification: {result}\nReason: {reason}\n")