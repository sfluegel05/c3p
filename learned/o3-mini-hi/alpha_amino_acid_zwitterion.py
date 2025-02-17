"""
Classifies: CHEBI:78608 alpha-amino-acid zwitterion
"""
#!/usr/bin/env python
"""
Classifies: alpha-amino-acid zwitterion
Definition: An amino acid‐zwitterion obtained by transfer of a proton from the carboxyl 
to the amino group of any alpha‐amino acid; major species at pH 7.3.
This version uses relaxed criteria to identify an alpha carbon that is bound to one 
protonated nitrogen and one deprotonated carboxylate group.
"""

from rdkit import Chem

def is_deprotonated_carboxylate(carbon_atom):
    """
    Checks if a given carbon atom represents a deprotonated carboxylate center (i.e. in -C(=O)[O-]).
    To be a deprotonated carboxylate carbon, we require:
      • The atom is carbon and sp2‐hybridized.
      • It is bonded to exactly two oxygen atoms.
      • One oxygen is double‐bonded with formal charge 0.
      • The other oxygen is single‐bonded with formal charge -1.
    We relax the requirement on the total heavy-atom count (to accommodate cyclic cases).
    """
    if carbon_atom.GetAtomicNum() != 6:
        return False
    if carbon_atom.GetHybridization() != Chem.rdchem.HybridizationType.SP2:
        return False
    # Get oxygen neighbors:
    oxy_neighbors = [nbr for nbr in carbon_atom.GetNeighbors() if nbr.GetAtomicNum() == 8]
    if len(oxy_neighbors) != 2:
        return False
    dbl_count = 0
    sgl_count = 0
    for bond in carbon_atom.GetBonds():
        nbr = bond.GetOtherAtom(carbon_atom)
        if nbr.GetAtomicNum() == 8:
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and nbr.GetFormalCharge() == 0:
                dbl_count += 1
            elif bond.GetBondType() == Chem.rdchem.BondType.SINGLE and nbr.GetFormalCharge() == -1:
                sgl_count += 1
    if dbl_count == 1 and sgl_count == 1:
        return True
    return False

def is_alpha_amino_acid_zwitterion(smiles: str):
    """
    Determines whether a molecule represents an alpha-amino-acid zwitterion.
    The heuristic is as follows:
      1) Parse the molecule and add explicit hydrogens.
      2) Iterate over candidate alpha carbons (aliphatic, sp3, with ≥1 hydrogen).
      3) For each candidate, check its heavy-atom neighbors for:
            a) Exactly one protonated nitrogen (atomic num 7, formal charge > 0, and at least one hydrogen).
            b) Exactly one deprotonated carboxylate carbon (recognized via is_deprotonated_carboxylate).
      4) If such a candidate is found, return True along with an explanation.
    If no such center is found, return False.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as an alpha-amino-acid zwitterion, False otherwise.
        str: Explanation for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure explicit hydrogens are added to properly count attached H atoms.
    mol = Chem.AddHs(mol)
    
    # Iterate over candidate alpha carbons.
    for atom in mol.GetAtoms():
        # Consider only carbon atoms.
        if atom.GetAtomicNum() != 6:
            continue
        # The alpha carbon should be sp3.
        if atom.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
            continue
        # Must have at least one bound hydrogen.
        if atom.GetTotalNumHs() < 1:
            continue
        
        # Check neighbors for one protonated nitrogen and one carboxylate carbon.
        protonated_nitrogen = None
        carboxylate_carbon = None
        
        # Evaluate heavy-atom neighbors.
        heavy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
        for nbr in heavy_neighbors:
            # Check for protonated nitrogen.
            if nbr.GetAtomicNum() == 7 and nbr.GetFormalCharge() > 0:
                # Verify that this nitrogen has at least one hydrogen.
                if nbr.GetTotalNumHs() >= 1:
                    if protonated_nitrogen is None:
                        protonated_nitrogen = nbr
                    else:
                        # If more than one protonated nitrogen neighbor, skip this candidate.
                        protonated_nitrogen = None
                        break
            # Check for deprotonated carboxylate carbon.
            if nbr.GetAtomicNum() == 6:
                if is_deprotonated_carboxylate(nbr):
                    if carboxylate_carbon is None:
                        carboxylate_carbon = nbr
                    else:
                        # More than one carboxylate group attached; not typical for a simple amino acid.
                        carboxylate_carbon = None
                        break
        
        # If both a unique protonated nitrogen and unique carboxylate carbon are found, we classify it as an alpha-amino-acid zwitterion.
        if protonated_nitrogen is not None and carboxylate_carbon is not None:
            return True, ("Found candidate alpha carbon (atom idx {}) attached to protonated nitrogen (atom idx {}) "
                          "and deprotonated carboxylate (atom idx {})."
                          .format(atom.GetIdx(), protonated_nitrogen.GetIdx(), carboxylate_carbon.GetIdx()))
    return False, "No alpha carbon with both connected deprotonated carboxylate and protonated amine found"

# Example usage:
if __name__ == "__main__":
    # These are a few examples taken from the provided list. In actual testing, one would add all examples.
    test_examples = [
        "[C@H]1(CCC[NH+]1C)C([O-])=O",               # N-methylproline zwitterion
        "[C@@H]([C@@H](C([O-])=O)[NH3+])(O)C#C",       # L-beta-ethynylserine zwitterion
        "CC[C@H](C)[C@H]([NH3+])C([O-])=O",            # L-isoleucine zwitterion
        "[NH3+][C@@H](CO)C([O-])=O",                   # L-serine zwitterion
    ]
    for smi in test_examples:
        res, msg = is_alpha_amino_acid_zwitterion(smi)
        print("SMILES:", smi)
        print("Result:", res)
        print("Reason:", msg)
        print("---------------------------")