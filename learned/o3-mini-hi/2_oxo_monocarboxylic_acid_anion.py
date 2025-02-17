"""
Classifies: CHEBI:35179 2-oxo monocarboxylic acid anion
"""
"""
Classifies: 2-oxo monocarboxylic acid anion
Definition: An oxo monocarboxylic acid anion in which the oxo (carbonyl) group is located at the 2‑position.
That is, if you number the chain from the terminal carboxylate (–C(=O)[O–], taken as position 1),
then the α‑position (position 2) should carry a carbonyl function (either directly on that carbon or via one of its substituents).
The algorithm first eliminates molecules containing common metal atoms to avoid false positives.
Then it searches for deprotonated carboxylate groups ("[CX3](=O)[O-]") that are terminal (attached to exactly one carbon).
Finally, it inspects the neighboring (α‑) carbon for the presence of a neutral carbonyl function.
"""

from rdkit import Chem

def is_2_oxo_monocarboxylic_acid_anion(smiles: str):
    """
    Determines whether a molecule is a 2-oxo monocarboxylic acid anion, meaning that
    it contains a terminal carboxylate group ([C](=O)[O-]) whose single carbon neighbor
    (the α‑carbon) either is itself a carbonyl carbon or has a substituent that is a carbonyl group.
    
    Args:
        smiles (str): SMILES string representing the molecule.
        
    Returns:
        bool: True if the molecule meets the definition.
        str: Explanation for the decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Reject molecules that contain metal atoms (commonly seen in false positives).
    # Here we exclude common metals that have been shown to give false positives.
    metals = {"Li", "Na", "K", "Mg", "Al", "Ca", "Fe", "Cu", "Zn"}
    for atom in mol.GetAtoms():
        if atom.GetSymbol() in metals:
            return False, "Molecule contains metal atoms, likely a salt or coordination complex"

    # Find candidate deprotonated (anionic) carboxylate groups using a SMARTS pattern.
    # The pattern "[CX3](=O)[O-]" matches a carboxylate carbon with one double-bonded oxygen and one negatively charged oxygen.
    carboxylate_smarts = "[CX3](=O)[O-]"
    carbox_mol = Chem.MolFromSmarts(carboxylate_smarts)
    carbox_matches = mol.GetSubstructMatches(carbox_mol)
    
    if not carbox_matches:
        return False, "No carboxylate group found"
    
    # Loop over each match and check if it forms a valid terminal carboxylate motif.
    for match in carbox_matches:
        # match[0] is the carboxylate carbon (the [C] from the pattern).
        carbox_c_idx = match[0]
        carbox_atom = mol.GetAtomWithIdx(carbox_c_idx)
        
        # Retrieve the neighboring carbons of the carboxylate carbon.
        carbon_neighbors = [nbr for nbr in carbox_atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
        
        # A terminal carboxylate should have exactly one carbon neighbor (the α‑carbon).
        if len(carbon_neighbors) != 1:
            continue  # Not terminal; try the next match.
        alpha_atom = carbon_neighbors[0]
        
        # Check if the α‑carbon itself bears a carbonyl; examine bonds from alpha_atom.
        has_carbonyl = False
        for bond in alpha_atom.GetBonds():
            # Get the bond type as a float; RDKit represents a double bond as 2.
            if bond.GetBondTypeAsDouble() == 2:
                # Identify the other atom in the double bond.
                other = bond.GetOtherAtom(alpha_atom)
                if other.GetAtomicNum() == 8 and other.GetFormalCharge() == 0:
                    has_carbonyl = True
                    break
        
        # If the α‑carbon is not directly a carbonyl, check whether any of its substituents (other than the terminal carboxylate carbon)
        # carries a carbonyl group.
        if not has_carbonyl:
            for nbr in alpha_atom.GetNeighbors():
                # Exclude the terminal carboxylate itself.
                if nbr.GetIdx() == carbox_atom.GetIdx():
                    continue
                if nbr.GetAtomicNum() != 6:
                    continue
                # For each neighboring carbon, inspect its bonds for a double-bonded oxygen.
                for bond in nbr.GetBonds():
                    if bond.GetBondTypeAsDouble() == 2:
                        oxy = bond.GetOtherAtom(nbr)
                        if oxy.GetAtomicNum() == 8 and oxy.GetFormalCharge() == 0:
                            has_carbonyl = True
                            break
                if has_carbonyl:
                    break
        
        if has_carbonyl:
            reason = ("Contains a terminal carboxylate group with an α‑carbon bearing a keto functionality "
                      "(directly as a carbonyl or via a substituent), matching the 2-oxo monocarboxylic acid anion motif.")
            return True, reason
    # If none of the candidate carboxylate groups qualifies, return false.
    return False, ("Does not contain a terminal carboxylate group with a 2‑position oxo (carbonyl) function. "
                   "Either the carboxylate is not terminal or the α‑carbon lacks a keto substituent.")

# Example test (uncomment to try a few molecules):
# test_smiles = [
#     "C\\C=C/CC(=O)C([O-])=O",  # cis-2-oxohex-4-enoate (expected True)
#     "[O-]C(C(CS(=O)[O-])=O)=O",  # 3-sulfinatopyruvate(2-) (expected True)
#     "CC(=O)NCCCCC([O-])=O",  # 5-acetamidopentanoate (expected False: missed motif)
#     "C(=O)(CC([O-])=O)C.[O-]C(CC(=O)C)=O"  # example with metal salt fragments - should be rejected
# ]
#
# for smi in test_smiles:
#     res, msg = is_2_oxo_monocarboxylic_acid_anion(smi)
#     print(f"SMILES: {smi}\nResult: {res}\nReason: {msg}\n{'-'*60}")