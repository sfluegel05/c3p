"""
Classifies: CHEBI:78608 alpha-amino-acid zwitterion
"""
#!/usr/bin/env python
"""
Classifies: alpha-amino-acid zwitterion
Definition: An amino acid‐zwitterion obtained by transfer of a proton from the carboxy to the amino group
of any alpha‐amino acid; major species at pH 7.3.
This version first defines a helper to strictly recognize a deprotonated carboxylate group.
Then it iterates over candidate alpha carbons (aliphatic, sp3, and with at least one hydrogen)
and requires that the candidate is directly attached to exactly one carboxylate group (as defined)
and exactly one nitrogen atom that is protonated (formal charge > 0 and with at least one hydrogen).
If such a center is found, we report success.
Note: This heuristic (being rule‐based) is not perfect.
"""

from rdkit import Chem

def is_deprotonated_carboxylate(carbon_atom):
    """
    Checks whether a given carbon atom is a deprotonated carboxylate center.
    In our heuristic, a proper –COO– carbon should:
      • be sp2,
      • have exactly two oxygen neighbors,
      • one oxygen should be connected by a double bond with no formal charge 
        and the other by a single bond with a formal charge of –1,
      • and the only heavy-atom neighbor (other than those two oxygens) should be the candidate alpha carbon.
    """
    if carbon_atom.GetAtomicNum() != 6:
        return False
    # Ideally the carboxyl carbon is sp2
    if carbon_atom.GetHybridization() != Chem.rdchem.HybridizationType.SP2:
        return False
    # Count oxygen neighbors
    oxy_neighbors = [nbr for nbr in carbon_atom.GetNeighbors() if nbr.GetAtomicNum() == 8]
    if len(oxy_neighbors) != 2:
        return False
    dbl_count = 0
    sgl_count = 0
    for bond in carbon_atom.GetBonds():
        nbr = bond.GetOtherAtom(carbon_atom)
        if nbr.GetAtomicNum() == 8:
            # Check bond type and formal charge on oxygen
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and nbr.GetFormalCharge() == 0:
                dbl_count += 1
            elif bond.GetBondType() == Chem.rdchem.BondType.SINGLE and nbr.GetFormalCharge() == -1:
                sgl_count += 1
    if dbl_count == 1 and sgl_count == 1:
        # Also require that besides the two oxygens, there is exactly one other heavy-atom neighbor.
        heavy_neighbors = [nbr for nbr in carbon_atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
        if len(heavy_neighbors) == 3:  
            # (alpha carbon + the 2 oxygens)
            return True
    return False

def is_alpha_amino_acid_zwitterion(smiles: str):
    """
    Determines if a molecule is an alpha-amino-acid zwitterion based on its SMILES string.
    The method looks for a carbon (the candidate alpha carbon) that is:
       - aliphatic, sp3,
       - has at least one bound hydrogen,
       - and is directly bonded to exactly one deprotonated carboxylate (–C(=O)[O–]) 
         and exactly one protonated nitrogen (e.g. [NH3+] or, in cyclic cases, [NH2+]).
    This heuristic is intended to catch common representations (including those with chiral tags)
    while avoiding many false positives from peptides or related compounds.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if classified as an alpha-amino-acid zwitterion, False otherwise.
        str: Explanation string as the reason for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Pre-calculate explicit hydrogens:
    mol = Chem.AddHs(mol)
    
    # We will iterate over candidate alpha carbons:
    for atom in mol.GetAtoms():
        # Look only at carbon atoms
        if atom.GetAtomicNum() != 6:
            continue
        # Require sp3 hybridization; note: sometimes chiral centers are marked but hybridization remains SP3
        if atom.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
            continue
        # Must have at least one hydrogen attached (explicit or implicit)
        if atom.GetTotalNumHs() < 1:
            continue
        # Now examine the heavy-atom neighbors (ignoring hydrogens)
        heavy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
        # For an alpha carbon in an amino acid we expect it to be attached to:
        #   (a) one protonated nitrogen, and (b) one carboxylate carbon (that is part of a -C(=O)[O-] group).
        protonated_N = None
        carboxylate_C = None
        
        for nbr in heavy_neighbors:
            # Check for a protonated amine neighbor
            if nbr.GetAtomicNum() == 7 and nbr.GetFormalCharge() > 0:
                # Check that this nitrogen has at least one hydrogen (to exclude quaternary-like ammoniums)
                if nbr.GetTotalNumHs() >= 1:
                    # If more than one such nitrogen shows up, we do not want to be ambiguous.
                    if protonated_N is None:
                        protonated_N = nbr
                    else:
                        protonated_N = None
                        break
            # Check for deprotonated carboxylate neighbor:
            if nbr.GetAtomicNum() == 6:
                # To be part of a deprotonated carboxylate, the neighbor carbon should be sp2
                # and must pass our strict check.
                if is_deprotonated_carboxylate(nbr):
                    if carboxylate_C is None:
                        carboxylate_C = nbr
                    else:
                        # More than one carboxylate attached? Not what we expect for a single amino acid.
                        carboxylate_C = None
                        break
        # If we found exactly one protonated nitrogen and one carboxylate, decide yes.
        if protonated_N is not None and carboxylate_C is not None:
            return True, ("Found matching alpha carbon (atom idx {}) attached to deprotonated carboxylate (atom idx {}) "
                          "and protonated amine (atom idx {}).".format(atom.GetIdx(), carboxylate_C.GetIdx(), protonated_N.GetIdx()))
    
    return False, "No alpha carbon with both connected deprotonated carboxylate and protonated amine found"

# Example usage (you may comment these out when using this module as an import):
if __name__ == "__main__":
    # Test with some examples:
    test_smiles_list = [
        "[C@H]1(CCC[NH+]1C)C([O-])=O",   # N-methylproline zwitterion: expected True
        "[C@@H]([C@@H](C([O-])=O)[NH3+])(O)C#C",   # L-beta-ethynylserine zwitterion: expected True
        "C[NH2+]CC([O-])=O",  # sarcosine zwitterion: expected False under improved criteria
        "[O-]C(=O)C([NH3+])CC=1N(C=NC1)C",  # 3-methylhistidine zwitterion: expected True
    ]
    for smi in test_smiles_list:
        res, reason = is_alpha_amino_acid_zwitterion(smi)
        print("SMILES:", smi)
        print("Result:", res)
        print("Reason:", reason)
        print("---------------------------")