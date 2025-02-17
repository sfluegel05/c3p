"""
Classifies: CHEBI:35179 2-oxo monocarboxylic acid anion
"""
"""
Classifies: 2-oxo monocarboxylic acid anion
Definition: An oxo monocarboxylic acid anion in which the oxo (carbonyl) group is located at the 2‑position.
That is, if you number the chain from the terminal carboxylate (–C(=O)[O–], taken as position 1),
then the α‑position (position 2) should carry a carbonyl function (either as the main chain atom or as a substituent).
This version tries to improve over a simple SMARTS match by (i) identifying a terminal carboxylate and (ii)
inspecting its immediate neighbor.
"""

from rdkit import Chem

def is_2_oxo_monocarboxylic_acid_anion(smiles: str):
    """
    Determines if a molecule is a 2-oxo monocarboxylic acid anion (α-keto acid anion)
    based on its SMILES string. In our definition, the molecule must have a terminal
    carboxylate group ([C](=O)[O-]) and the carbon next to that group (the α‐carbon)
    must either itself be a carbonyl carbon or have a substituent that is a carbonyl group.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is classified as a 2-oxo monocarboxylic acid anion, False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # First, find terminal carboxylate groups.
    # The SMARTS "[CX3](=O)[O-]" should match a deprotonated carboxylate.
    carboxylate_smarts = "[CX3](=O)[O-]"
    carbox_mol = Chem.MolFromSmarts(carboxylate_smarts)
    carbox_matches = mol.GetSubstructMatches(carbox_mol)
    
    if not carbox_matches:
        return False, "No carboxylate group found"
    
    # We now require that ONE of the carboxylate groups is terminal:
    # i.e. its carboxylate carbon is attached to exactly one carbon (the α-carbon).
    valid_motif_found = False
    reason_details = ""
    
    for match in carbox_matches:
        # In our pattern, match[0] is the carboxylate carbon.
        carbox_c = match[0]
        atom_c = mol.GetAtomWithIdx(carbox_c)
        
        # Find all neighbors of this carbon that are carbons.
        neighbor_carbons = [nbr for nbr in atom_c.GetNeighbors() if nbr.GetAtomicNum() == 6]
        
        # We require that the carboxylate carbon has exactly one carbon neighbor.
        if len(neighbor_carbons) != 1:
            continue  # Not terminal: skip
        
        # We designate the sole carbon neighbor as the α‑carbon.
        alpha_atom = neighbor_carbons[0]
        
        # Now, check for the keto function at the α‐position.
        # It might be that the α‐carbon itself is a carbonyl C (has a C=O bond with a non‐anionic oxygen)
        has_ketone = False
        
        # Check if the alpha_atom is itself carbonyl (has a double bond to an oxygen that is not carrying a negative charge)
        for bond in alpha_atom.GetBonds():
            if bond.GetBondTypeAsDouble() == 2:
                nbr = bond.GetOtherAtom(alpha_atom)
                # Accept if oxygen and not [O-] (note: in typical keto groups the O is neutral)
                if nbr.GetAtomicNum() == 8 and nbr.GetFormalCharge() == 0:
                    has_ketone = True
                    break
        
        # If not directly a carbonyl, check if any substituent (other than the terminal carboxylate)
        # attached to the α‐carbon is a carbonyl carbon.
        if not has_ketone:
            for nbr in alpha_atom.GetNeighbors():
                # Exclude the terminal carboxylate carbon.
                if nbr.GetIdx() == atom_c.GetIdx():
                    continue
                if nbr.GetAtomicNum() != 6:
                    continue
                # Check if this neighbor has a double-bonded oxygen.
                for bond in nbr.GetBonds():
                    if bond.GetBondTypeAsDouble() == 2:
                        oxy = bond.GetOtherAtom(nbr)
                        if oxy.GetAtomicNum() == 8 and oxy.GetFormalCharge() == 0:
                            has_ketone = True
                            break
                if has_ketone:
                    break
        
        if has_ketone:
            valid_motif_found = True
            reason_details = ("Contains a terminal carboxylate group with an α‑carbon "
                              "bearing a keto function (directly or via a substituent), matching "
                              "the 2-oxo monocarboxylic acid anion motif.")
            break  # We only need one valid motif.
    
    if valid_motif_found:
        return True, reason_details
    else:
        return False, ("Does not contain a terminal carboxylate group with a 2‑position oxo (carbonyl) "
                       "function. Either the carboxylate is not terminal "
                       "or the α‑carbon lacks a keto substituent.")


# Uncomment below lines to run a few tests:
# test_smiles_list = [
#     "C\\C=C/CC(=O)C([O-])=O",  # cis-2-oxohex-4-enoate (should be True)
#     "[O-]C(C(CS(=O)[O-])=O)=O",  # 3-sulfinatopyruvate(2-) (should be True)
#     "C(NC(=[NH2+])N(C)C)CCC(C(=O)[O-])=O",  # false positive example should return False
#     "CC(=O)C(C)(O)C([O-])=O",  # 2-acetyllactate (expected True)
#     "CC(=O)NCCCCC([O-])=O",   # 5-acetamidopentanoate (expected True)
#     "OC[C@@H](O)[C@H](O)[C@H](O)CC(=O)C([O-])=O"  # example false positive from prior approach
# ]
# 
# for smi in test_smiles_list:
#     result, reason = is_2_oxo_monocarboxylic_acid_anion(smi)
#     print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n{'-'*60}")