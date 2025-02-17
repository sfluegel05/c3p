"""
Classifies: CHEBI:50760 N-hydroxy-alpha-amino-acid
"""
"""
Classifies: N-hydroxy-alpha-amino-acid
Definition: Any amino acid in which at least one hydrogen attached to the amino group is replaced by a hydroxy group.
In other words, we look for a typical α–amino acid backbone (carboxyl group attached to an α–carbon)
where the amino substituent has an –OH group.
"""
from rdkit import Chem

def is_N_hydroxy_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-hydroxy-alpha-amino-acid based on its SMILES string.
    
    It does so by:
      1. Checking that a carboxyl group (C(=O)O) is present.
      2. Identifying the α–carbon as the non-oxygen neighbor of the carboxyl carbon.
      3. Verifying that the α–carbon is bonded to a nitrogen (the amino group)
         that carries at least one hydroxy (-OH) substituent.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as an N-hydroxy-alpha-amino-acid, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string to obtain a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Define a SMARTS pattern for a carboxyl group (C(=O)O)
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    if not carboxyl_matches:
        return False, "No carboxyl group found; not a typical amino acid backbone."
    
    # For each carboxyl group match, try to find the alpha carbon and then check for a hydroxy-substituted nitrogen.
    for match in carboxyl_matches:
        # In the SMARTS "C(=O)O", match[0] is the carboxyl carbon.
        carboxyl_carbon_idx = match[0]
        carboxyl_carbon = mol.GetAtomWithIdx(carboxyl_carbon_idx)
        
        # Look for the alpha carbon: it should be a neighbor of the carboxyl carbon that is not an oxygen.
        alpha_carbons = []
        for neighbor in carboxyl_carbon.GetNeighbors():
            if neighbor.GetSymbol() != "O":
                alpha_carbons.append(neighbor)
        
        if not alpha_carbons:
            continue  # try next match; this one does not lead to a potential α–carbon
        
        # For each candidate alpha carbon, look for an amino (nitrogen) neighbor (excluding the carboxyl carbon)
        for alpha_carbon in alpha_carbons:
            # Check neighbors of the alpha carbon (excluding the carboxyl carbon)
            for nbr in alpha_carbon.GetNeighbors():
                if nbr.GetIdx() == carboxyl_carbon_idx:
                    continue
                # We are interested in nitrogen atoms (the amino group)
                if nbr.GetSymbol() == "N":
                    amino_nitrogen = nbr
                    # Check if the nitrogen has at least one oxygen substituent attached via a single bond.
                    for oxy_neighbor in amino_nitrogen.GetNeighbors():
                        # We consider the O substituent only if it is not the alpha carbon itself.
                        if oxy_neighbor.GetSymbol() == "O" and oxy_neighbor.GetIdx() != alpha_carbon.GetIdx():
                            # Check that the bond between the nitrogen and the oxygen is a single bond
                            bond = mol.GetBondBetweenAtoms(amino_nitrogen.GetIdx(), oxy_neighbor.GetIdx())
                            if bond and bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                                return True, "Found α–carbon with amino nitrogen carrying an –OH substituent."
    # If we did not find any amino nitrogen with an –OH attached to an alpha carbon, return False.
    return False, "Molecule does not present an N-hydroxy amino group attached to an α–carbon in a carboxylate context."

# For testing purposes (optional)
if __name__ == "__main__":
    test_smiles = [
        "O=C(O)[C@@H](NO)CCCCSC",         # N-hydroxy-L-dihomomethionine
        "CSCCCCCCCCC(NO)C(O)=O",            # N-hydroxyhexahomomethionine
        "ONCC(O)=O",                       # N-hydroxyglycine
        "CC(C)[C@H](N(O)O)C(O)=O",          # N,N-dihydroxy-L-valine (should be True as well)
        "CC(C)[C@H](N)C(O)=O",              # L-valine without N-hydroxy (should be False)
    ]
    
    for smi in test_smiles:
        is_class, reason = is_N_hydroxy_alpha_amino_acid(smi)
        print(f"SMILES: {smi}\nClassification: {is_class}\nReason: {reason}\n")