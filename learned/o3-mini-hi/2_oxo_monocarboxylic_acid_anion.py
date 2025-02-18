"""
Classifies: CHEBI:35179 2-oxo monocarboxylic acid anion
"""
"""
Classifies: 2-oxo monocarboxylic acid anion
Definition: An oxo monocarboxylic acid anion in which the oxo (carbonyl) group is located at the 2‑position.
That is, if you number the chain from the terminal carboxylate (–C(=O)[O–], taken as position 1),
then the α‑position (position 2) must bear a carbonyl functionality (directly, via a double-bond, 
or via a substituent).
This version improves upon a previous attempt by ensuring that the matched carboxylate group is terminal 
(i.e. the carboxylate carbon is attached only to two oxygen atoms and exactly one carbon “α‑carbon”) 
and then checking that the α‑carbon carries a neutral carbonyl group.
"""

from rdkit import Chem

def is_2_oxo_monocarboxylic_acid_anion(smiles: str):
    """
    Determines whether a molecule is a 2-oxo monocarboxylic acid anion.
    This requires that the molecule contains a terminal carboxylate group ([C](=O)[O-]) where the carboxylate
    carbon is attached to exactly one carbon neighbor (the α‑carbon) and that this α‑carbon carries a carbonyl 
    functionality (i.e. has a double bond to oxygen with formal charge zero) in one of its substituents.
    
    Args:
        smiles (str): SMILES string representing the molecule.
    
    Returns:
        bool: True if the molecule satisfies the definition.
        str: Explanation for the decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Exclude molecules containing common metal atoms – these are often salts or complexes
    metals = {"Li", "Na", "K", "Mg", "Al", "Ca", "Fe", "Cu", "Zn"}
    for atom in mol.GetAtoms():
        if atom.GetSymbol() in metals:
            return False, "Molecule contains metal atoms, likely a salt or coordination complex"
    
    # Find candidate carboxylate groups using SMARTS.
    # This pattern matches a carboxylate carbon with one carbonyl oxygen and one negatively charged oxygen.
    carboxylate_smarts = "[CX3](=O)[O-]"
    carbox_mol = Chem.MolFromSmarts(carboxylate_smarts)
    carbox_matches = mol.GetSubstructMatches(carbox_mol)
    
    if not carbox_matches:
        return False, "No carboxylate group found"
    
    # Process each matched carboxylate candidate.
    for match in carbox_matches:
        # match[0] is the carboxylate carbon (the [C] in the pattern).
        carbox_c_idx = match[0]
        carbox_atom = mol.GetAtomWithIdx(carbox_c_idx)
        
        # Check that the carboxylate carbon has exactly three heavy-atom neighbors:
        # two oxygens (the carbonyl and the negatively charged oxygen) and one carbon.
        neighbors = [nbr for nbr in carbox_atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
        if len(neighbors) != 3:
            continue  # Not a terminal carboxylate; try next candidate
        
        # Identify the single carbon neighbor (the α‑carbon)
        carbon_neighbors = [nbr for nbr in neighbors if nbr.GetAtomicNum() == 6]
        if len(carbon_neighbors) != 1:
            continue  # More than one carbon neighbor disqualifies the candidate
        alpha_atom = carbon_neighbors[0]
        
        # Now check the α‑carbon for carbonyl functionality.
        # We need to ensure that aside from being bonded to the carboxylate carbon,
        # at least one bond emanating from the α‑carbon is a double bond to oxygen (i.e. a carbonyl)
        has_carbonyl = False
        
        # Look directly on the α‑carbon (excluding the bond back to the carboxylate group)
        for bond in alpha_atom.GetBonds():
            # Identify the bond partner.
            other_atom = bond.GetOtherAtom(alpha_atom)
            if other_atom.GetIdx() == carbox_atom.GetIdx():
                continue  # skip the carboxylate connection
            if bond.GetBondTypeAsDouble() == 2:
                # It’s a double bond; check if the partner atom is oxygen and is neutral.
                if other_atom.GetAtomicNum() == 8 and other_atom.GetFormalCharge() == 0:
                    has_carbonyl = True
                    break
        
        # Optionally, one might search one step further (i.e. a substituent on the α‑carbon that is itself a carbonyl)
        # but this simple direct check already greatly reduces false positives.
        
        if has_carbonyl:
            reason = ("Contains a terminal carboxylate group (C(=O)[O–]) attached to an α‑carbon that bears a carbonyl "
                      "functionality, consistent with the 2‑oxo monocarboxylic acid anion definition.")
            return True, reason

    # If no candidate matches the criteria, then it is either missing a terminal carboxylate or the required α‑carbon carbonyl.
    return False, ("Does not contain a terminal carboxylate group with an α‑carbon that bears a carbonyl functionality "
                   "at the 2‑position.")

# Example test (uncomment the lines below to run tests):
# test_smiles = [
#     "C\\C=C/CC(=O)C([O-])=O",  # cis-2-oxohex-4-enoate (expected True)
#     "[O-]C(C(CS(=O)[O-])=O)=O",  # 3-sulfinatopyruvate(2-) (expected True)
#     "CC(=O)NCCCCC([O-])=O",  # 5-acetamidopentanoate (expected False, missing α‑carbon carbonyl)
#     "OC(=O)C([O-])=O",  # glyoxylate (expected True: terminal carboxylate and α‑position is the carbonyl of glyoxylate itself)
# ]
# 
# for smi in test_smiles:
#     res, msg = is_2_oxo_monocarboxylic_acid_anion(smi)
#     print(f"SMILES: {smi}\nResult: {res}\nReason: {msg}\n{'-'*60}")