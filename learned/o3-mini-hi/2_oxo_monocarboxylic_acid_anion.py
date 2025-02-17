"""
Classifies: CHEBI:35179 2-oxo monocarboxylic acid anion
"""
"""
Classifies: 2-oxo monocarboxylic acid anion
Definition: An oxo monocarboxylic acid anion in which the oxo (carbonyl) group is located at the 2‑position.
That is, when numbering the chain from a terminal carboxylate group (–C(=O)[O–], where the carboxylate carbon is 
bonded to exactly two oxygens and exactly one carbon), the unique α‑carbon must bear a carbonyl functionality.
This functionality may be either on the α‑carbon itself (a direct C=O) or on a substituent (an acyl group bound 
to an sp³ α‑carbon). Only molecules with a single terminal carboxylate group are accepted.
"""

from rdkit import Chem

def is_2_oxo_monocarboxylic_acid_anion(smiles: str):
    """
    Determines whether a molecule is a 2-oxo monocarboxylic acid anion.
    The molecule must have one (and only one) terminal carboxylate group 
    ([CX3](=O)[O-]) where the carboxylate carbon is attached to exactly two oxygens 
    (one via a double bond, one with a negative charge) and exactly one carbon neighbor 
    (the α‑carbon). The α‑carbon must either directly have a carbonyl (a double bond to oxygen, formal charge 0)
    or have a substituent (an acyl group) in which the attached carbon carries a carbonyl.
    
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
    
    # Exclude molecules that contain common metals (often salts/complexes).
    metals = {"Li", "Na", "K", "Mg", "Al", "Ca", "Fe", "Cu", "Zn"}
    for atom in mol.GetAtoms():
        if atom.GetSymbol() in metals:
            return False, "Molecule contains metal atoms; likely not a valid organic monocarboxylate"
    
    # Find candidate terminal carboxylate groups using SMARTS.
    # This SMARTS matches a carboxylate carbon (bonded to one double-bonded oxygen and one negatively charged oxygen).
    carboxylate_smarts = "[CX3](=O)[O-]"
    carbox_mol = Chem.MolFromSmarts(carboxylate_smarts)
    carbox_matches = mol.GetSubstructMatches(carbox_mol)
    
    # Keep only those candidates that are terminal: carboxyl carbon must have exactly 3 heavy neighbors 
    # (2 oxygens and 1 carbon).
    terminal_candidates = []
    for match in carbox_matches:
        # match[0] corresponds to the carboxylate carbon.
        c_idx = match[0]
        carbox_c = mol.GetAtomWithIdx(c_idx)
        heavy_neighbors = [nbr for nbr in carbox_c.GetNeighbors() if nbr.GetAtomicNum() > 1]
        # Count oxygen atoms.
        oxy_neighbors = [nbr for nbr in heavy_neighbors if nbr.GetAtomicNum() == 8]
        carbon_neighbors = [nbr for nbr in heavy_neighbors if nbr.GetAtomicNum() == 6]
        # Terminal means: exactly 2 oxygen neighbors and exactly 1 carbon neighbor.
        if len(oxy_neighbors) == 2 and len(carbon_neighbors) == 1:
            terminal_candidates.append((c_idx, carbon_neighbors[0].GetIdx()))
    
    # We require exactly ONE terminal carboxylate group.
    if len(terminal_candidates) == 0:
        return False, "No terminal carboxylate group found"
    if len(terminal_candidates) > 1:
        return False, "Multiple terminal carboxylate groups found; not a monocarboxylic acid anion"
    
    # For the unique candidate, get the α‑carbon index.
    carbox_c_idx, alpha_idx = terminal_candidates[0]
    alpha_atom = mol.GetAtomWithIdx(alpha_idx)
    carbox_atom = mol.GetAtomWithIdx(carbox_c_idx)
    
    # Now check the α‑carbon for a carbonyl functionality.
    # Option (A): Direct check on α‑carbon for any double bond to oxygen (with formal charge 0)
    direct_carbonyl = False
    for bond in alpha_atom.GetBonds():
        # Skip the bond going back to the carboxylate carbon.
        if bond.GetOtherAtom(alpha_atom).GetIdx() == carbox_c_idx:
            continue
        if bond.GetBondTypeAsDouble() == 2:
            other = bond.GetOtherAtom(alpha_atom)
            if other.GetAtomicNum() == 8 and other.GetFormalCharge() == 0:
                direct_carbonyl = True
                break
    if direct_carbonyl:
        reason = ("Contains a unique terminal carboxylate group (C(=O)[O–]) attached to an α‑carbon that "
                  "bears a direct carbonyl (C=O), consistent with the 2‑oxo monocarboxylic acid anion definition.")
        return True, reason
    
    # Option (B): The α‑carbon might have an acyl substituent.
    # For each neighbor (other than the carboxylate group), check if it is a carbon that itself is carbonyl.
    acyl_carbonyl = False
    for nbr in alpha_atom.GetNeighbors():
        if nbr.GetIdx() == carbox_c_idx:
            continue
        # We require the neighbor is carbon.
        if nbr.GetAtomicNum() != 6:
            continue
        # Look for a double bond from this neighbor to an oxygen (with formal charge 0).
        for bond in nbr.GetBonds():
            # Allow the acyl carbon to be double-bonded to oxygen.
            if bond.GetBondTypeAsDouble() == 2:
                other = bond.GetOtherAtom(nbr)
                if other.GetAtomicNum() == 8 and other.GetFormalCharge() == 0:
                    acyl_carbonyl = True
                    break
        if acyl_carbonyl:
            break

    if acyl_carbonyl:
        reason = ("Contains a unique terminal carboxylate group (C(=O)[O–]) attached to an α‑carbon that bears an acyl substituent "
                  "with a carbonyl (C=O), consistent with the 2‑oxo monocarboxylic acid anion definition.")
        return True, reason

    # Neither a direct carbonyl on the α‑carbon nor an acyl substituent was found.
    return False, ("Does not contain a terminal carboxylate group with an α‑carbon that bears a carbonyl functionality "
                   "at the 2‑position.")
    
# Example tests (uncomment to run)
# test_smiles = [
#     "C\\C=C/CC(=O)C([O-])=O",  # cis-2-oxohex-4-enoate (expected True)
#     "[O-]C(C(CS(=O)[O-])=O)=O",  # 3-sulfinatopyruvate(2-) (expected True)
#     "OC1CCC(CC(=O)C([O-])=O)C=C1",  # tetrahydro-4-hydroxyphenylpyruvate (expected True)
#     "CC(=O)C(C)(O)C([O-])=O",  # 2-acetyllactate (expected True)
#     "CC(=O)NCCCCC([O-])=O",  # 5-acetamidopentanoate (expected False)
#     "[H]C(=O)[C@H](O)[C@@H](O)CC(=O)C([O-])=O",  # 5-dehydro-4-deoxy-D-glucuronate (expected False)
# ]
#
# for smi in test_smiles:
#     res, msg = is_2_oxo_monocarboxylic_acid_anion(smi)
#     print(f"SMILES: {smi}\nResult: {res}\nReason: {msg}\n{'-'*60}")