"""
Classifies: CHEBI:35179 2-oxo monocarboxylic acid anion
"""
"""
Classifies: 2-oxo monocarboxylic acid anion
Definition: An oxo monocarboxylic acid anion in which the oxo (carbonyl) group is located at the 2‑position.
That is, when numbering the chain from a terminal carboxylate group (–C(=O)[O–], where the carboxylate carbon is 
bonded to exactly two oxygens and exactly one carbon), the unique α‑carbon must:
  (A) either directly bear a double bond to oxygen (a carbonyl), or 
  (B) have exactly one non-carboxylate carbon substituent that itself is carbonyl.
Only molecules with exactly one terminal carboxylate group are accepted.
"""

from rdkit import Chem

def is_2_oxo_monocarboxylic_acid_anion(smiles: str):
    """
    Determines whether a molecule is a 2‑oxo monocarboxylic acid anion.
    It requires a unique terminal carboxylate group ([CX3](=O)[O-]) in which the carboxylate
    carbon is attached to exactly two oxygens and one carbon. That unique carbon neighbor (the α‑carbon)
    must bear a carbonyl at the 2‑position (either directly or via an acyl substituent).
    
    Args:
        smiles (str): SMILES string representing the molecule.
    
    Returns:
        bool: True if the molecule satisfies the definition.
        str: Explanation for the decision.
    """
    # Parse SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # (Optional) Exclude compounds with common metals.
    metals = {"Li", "Na", "K", "Mg", "Al", "Ca", "Fe", "Cu", "Zn"}
    for atom in mol.GetAtoms():
        if atom.GetSymbol() in metals:
            return False, "Molecule contains metal atoms; likely not a relevant organic monocarboxylate"
    
    # Find candidate terminal carboxylate groups using a SMARTS pattern.
    # This pattern matches a carboxylate carbon that is double-bonded to an oxygen and single-bonded to an [O-].
    carboxylate_smarts = "[CX3](=O)[O-]"
    carbox_mol = Chem.MolFromSmarts(carboxylate_smarts)
    carbox_matches = mol.GetSubstructMatches(carbox_mol)
    
    # Only consider candidates in which the carboxylate C has exactly two oxygen neighbors and exactly one carbon neighbor.
    terminal_candidates = []
    for match in carbox_matches:
        c_idx = match[0]  # index of the carboxylate carbon
        carbox_c = mol.GetAtomWithIdx(c_idx)
        heavy_neighbors = [nbr for nbr in carbox_c.GetNeighbors() if nbr.GetAtomicNum() > 1]
        oxy_neighbors = [nbr for nbr in heavy_neighbors if nbr.GetAtomicNum() == 8]
        carbon_neighbors = [nbr for nbr in heavy_neighbors if nbr.GetAtomicNum() == 6]
        if len(oxy_neighbors) == 2 and len(carbon_neighbors) == 1:
            terminal_candidates.append((c_idx, carbon_neighbors[0].GetIdx()))
    
    if len(terminal_candidates) == 0:
        return False, "No terminal carboxylate group found"
    if len(terminal_candidates) > 1:
        return False, "Multiple terminal carboxylate groups found; not a monocarboxylic acid anion"
    
    # We now have a unique terminal carboxylate candidate.
    carbox_c_idx, alpha_idx = terminal_candidates[0]
    alpha_atom = mol.GetAtomWithIdx(alpha_idx)
    
    # Option (A): Check if the α‑carbon itself bears a carbonyl.
    # That is, one of the bonds from the α‑carbon (except the one leading back to the carboxylate)
    # is a double bond to an oxygen (with formal charge 0). This means the α‑carbon itself is oxo.
    for bond in alpha_atom.GetBonds():
        other = bond.GetOtherAtom(alpha_atom)
        # Skip the bond back to the terminal carboxylate.
        if other.GetIdx() == carbox_c_idx:
            continue
        if bond.GetBondTypeAsDouble() == 2 and other.GetAtomicNum() == 8 and other.GetFormalCharge() == 0:
            reason = ("Contains a unique terminal carboxylate group (C(=O)[O–]) attached to an α‑carbon which "
                      "directly bears a carbonyl (C=O), consistent with the 2‑oxo monocarboxylic acid anion definition.")
            return True, reason
    
    # Option (B): Check if the α‑carbon bears an acyl substituent.
    # Here we require that among the α‑carbon’s neighbors (other than the terminal carboxylate),
    # there is exactly one carbon (acyl carbon) that is itself carbonyl (i.e. double‐bonded to oxygen with FC 0).
    acyl_candidates = []
    for nbr in alpha_atom.GetNeighbors():
        if nbr.GetIdx() == carbox_c_idx:
            continue
        # Only consider carbon neighbors.
        if nbr.GetAtomicNum() != 6:
            continue
        # Check if this neighbor has a double bond to oxygen (avoiding the bond going to α‑carbon).
        for bond in nbr.GetBonds():
            other = bond.GetOtherAtom(nbr)
            if other.GetIdx() == alpha_atom.GetIdx():
                continue
            if bond.GetBondTypeAsDouble() == 2 and other.GetAtomicNum() == 8 and other.GetFormalCharge() == 0:
                acyl_candidates.append(nbr.GetIdx())
                # Break out once we have found that this neighbor is acyl–like.
                break
    if len(acyl_candidates) == 1:
        reason = ("Contains a unique terminal carboxylate group (C(=O)[O–]) attached to an α‑carbon that bears an acyl substituent "
                  "with a carbonyl (C=O), consistent with the 2‑oxo monocarboxylic acid anion definition.")
        return True, reason

    return False, ("Does not contain a unique terminal carboxylate group with an α‑carbon that bears a carbonyl functionality "
                   "at the 2‑position.")

# Example tests (uncomment to run)
# test_smiles = [
#     "C\\C=C/CC(=O)C([O-])=O",  # cis-2-oxohex-4-enoate (expected True)
#     "[O-]C(C(CS(=O)[O-])=O)=O",  # 3-sulfinatopyruvate(2-) (expected True)
#     "OC1CCC(CC(=O)C([O-])=O)C=C1",  # tetrahydro-4-hydroxyphenylpyruvate (expected True)
#     "C(C([C@H](O)C)=O)(=O)[O-]",  # (R)-3-hydroxy-2-oxobutanoate (expected True)
#     "[O-]C(=O)C(=O)c1ccccc1",  # phenylglyoxylate (expected True)
#     "CC(=O)NCCCCC([O-])=O",  # 5-acetamidopentanoate (expected False per definition)
# ]
#
# for smi in test_smiles:
#     res, msg = is_2_oxo_monocarboxylic_acid_anion(smi)
#     print(f"SMILES: {smi}\nResult: {res}\nReason: {msg}\n{'-'*60}")