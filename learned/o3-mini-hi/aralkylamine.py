"""
Classifies: CHEBI:18000 aralkylamine
"""
"""
Classifies: aralkylamine
Definition: An alkylamine in which at least one substituent (attached through one or more bonds)
is an aromatic group. In our operational definition we require that the amine nitrogen be non‐aromatic
and that some aromatic carbon (assumed to be the marker for an aromatic ring) is found within a short bond-distance
(1 to 3 bonds) from that nitrogen.
"""

from rdkit import Chem
import numpy as np

def is_aralkylamine(smiles: str):
    """
    Determines if a molecule is an aralkylamine based on its SMILES string.
    
    Here our approach is:
      1. Parse the SMILES into an RDKit molecule.
      2. Identify all nitrogen atoms that are non‐aromatic 
         (i.e. considered “aliphatic”).
      3. Identify all aromatic carbon atoms (as markers for aromatic groups).
      4. Compute the bond-count distance matrix.
      5. For each non‐aromatic nitrogen, if any aromatic carbon is found 
         within a short path length (1–3 bonds) then we conclude that 
         the nitrogen is substituted (directly or via a short alkyl linker)
         by an aromatic group and classify the molecule as an aralkylamine.
      
    Note: This “distance‐based” rule is one way to balance the desire to catch short alkyl linkers 
    (benzylamine, phenethylamine, 2‐anilinoethanol) while not classifying very long “remote” aromatic groups 
    (as may occur in peptides) as aralkylamines.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as an aralkylamine, False otherwise.
        str: An explanation for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get the distance matrix (bond counts) for the molecule.
    dmat = Chem.GetDistanceMatrix(mol)
    
    # Identify indices for non-aromatic nitrogen atoms.
    # (Note: In aromatic amines the nitrogen is marked as aromatic so we skip those.)
    na_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7 and not atom.GetIsAromatic()]
    
    # Identify indices for aromatic carbons as markers for aromatic rings.
    ar_c_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetIsAromatic()]
    
    # If no non-aromatic nitrogens or no aromatic carbons, we cannot have an aralkylamine.
    if not na_indices or not ar_c_indices:
        return False, "No non‐aromatic amine or aromatic ring found in the molecule."
    
    # Define maximum allowed bond separation.
    # This value can be tuned: benzylamine has distance=2, phenethylamine has distance=3;
    # we consider that up to 3 bonds is a plausible limit for an alkyl chain substituent.
    max_distance = 3

    # Loop over each non-aromatic nitrogen.
    for n_idx in na_indices:
        for ar_idx in ar_c_indices:
            # dmat entry gives minimum path length (in bonds) between two atoms.
            dist = dmat[n_idx, ar_idx]
            # We require that the distance is at least 1 but no more than max_distance.
            if 1 <= dist <= max_distance:
                return True, ("Found a non‐aromatic amine nitrogen with an aromatic group "
                              f"within {dist} bond(s). Molecule classified as aralkylamine.")
    
    return False, ("No aralkylamine substructure found: no non‐aromatic amine is connected "
                   "to an aromatic group within 3 bonds.")

# Example test cases (these may be removed or commented out when used as a module)
if __name__ == "__main__":
    # Some examples – note that in our working definition:
    # * benzylamine (NCc1ccccc1) gives a path: N-CH2 (distance 1) then to aromatic C (distance 2).
    # * phenethylamine (NCCc1ccccc1) gives a path of 3 bonds.
    # * aniline (c1ccc(N)cc1) is not considered because the N is aromatic.
    # * 2-Anilinoethanol (OCCNC1=CC=CC=C1) — here the nitrogen is non-aromatic even though one substituent is directly aromatic,
    #   and will be classified (distance = 1).
    test_smiles = [
        ("NCc1ccccc1", "benzylamine"),
        ("NCCc1ccccc1", "phenethylamine"),
        ("c1ccc(N)cc1", "aniline (should not classify, N aromatic)"),
        ("OCCNC1=CC=CC=C1", "2-Anilinoethanol")
    ]
    
    for sm, name in test_smiles:
        result, reason = is_aralkylamine(sm)
        print(f"SMILES: {sm}\nMolecule: {name}\nClassification: {result}\nReason: {reason}\n")