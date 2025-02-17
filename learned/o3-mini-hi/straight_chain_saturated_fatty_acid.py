"""
Classifies: CHEBI:39418 straight-chain saturated fatty acid
"""
#!/usr/bin/env python
"""
Classifies: Straight-chain saturated fatty acid
Definition: Any saturated fatty acid lacking a side-chain.
Rules implemented here:
  (1) Molecule must contain exactly one carboxyl group (allowing both COOH and COO– forms).
  (2) When considering only the carbon atoms connected via single bonds, the carbons must form 
      a single linear (unbranched) chain (i.e. exactly two endpoints and all others of degree 2). 
  (3) The carboxyl carbon must be an endpoint.
  (4) For each carbon of the chain, substituents outside the chain are allowed only if:
        - For the carboxyl carbon: only the two oxygens that make up the carboxyl.
        - For the non-carboxyl terminal carbon: at most one –OH (i.e. one oxygen single-bonded) is allowed.
        - For an internal chain carbon: only a carbonyl oxygen (bonded by a double bond) is tolerated.
  (5) No unsaturation is allowed in the C–C main chain.
Any violation leads to a rejection.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_straight_chain_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a straight-chain saturated fatty acid.
    The criteria are:
      - Exactly one carboxyl group (COOH or COO–) is present.
      - All carbons (atomic number 6) connected by single bonds form one linear chain.
      - The carboxyl carbon must be at one end of the chain.
      - Aside from the carboxyl group, extra substituents are not allowed.
         (A terminal (non-carboxyl) carbon may bear at most one -OH group;
          an internal carbon may have a double-bonded oxygen (keto) group but no single-bonded extra substituents.)
      - All C–C bonds are single (no unsaturation in the main chain).
      
    Args:
         smiles (str): SMILES string of the molecule.
         
    Returns:
         (bool, str): True with an explanation if the molecule is classified as a straight‐chain saturated fatty acid;
                      False with a reason otherwise.
    """
    # Parse SMILES and add explicit hydrogens (to catch substituents such as –OH).
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)
    
    # (D) Check that there are no unsaturated C–C bonds.
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
            if bond.GetBondType() != Chem.BondType.SINGLE:
                return False, "Found an unsaturated carbon–carbon bond"
    
    # (A) Identify the carboxyl group.
    # SMARTS for carboxyl: a carbon with a double-bonded oxygen and a single-bonded oxygen (which may be protonated or be an anion).
    carboxyl_smarts = "C(=O)[O;H1,-1]"  # works for both COOH and COO-
    carboxyl_pattern = Chem.MolFromSmarts(carboxyl_smarts)
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    if len(carboxyl_matches) != 1:
        return False, f"Found {len(carboxyl_matches)} carboxyl groups, expected exactly 1"
    # Use the first match, where index 0 is the carboxyl carbon.
    carboxyl_match = carboxyl_matches[0]
    carboxyl_carbon_idx = carboxyl_match[0]
    allowed_for_carboxyl = set(carboxyl_match[1:])  # these oxygens are allowed on the carboxyl carbon

    # (B) Build a connectivity graph for all carbon atoms using only single bonds.
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    carbon_ids = [atom.GetIdx() for atom in carbon_atoms]
    # Create a dictionary to hold neighbors (only considering C–C single bonds).
    cgraph = {idx: [] for idx in carbon_ids}
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.SINGLE:
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
                i1 = a1.GetIdx()
                i2 = a2.GetIdx()
                if i1 in cgraph and i2 in cgraph:
                    cgraph[i1].append(i2)
                    cgraph[i2].append(i1)
    # The main carbon chain should include all carbon atoms.
    if set(cgraph.keys()) != set(carbon_ids):
        return False, "Not all carbons are in the connectivity graph"

    # Check that the carbon connectivity graph forms a linear chain (a simple path)
    # In a simple chain, there should be exactly two atoms with degree 1 and all others with degree 2.
    endpoints = [node for node, nbrs in cgraph.items() if len(nbrs)==1]
    if len(endpoints) != 2:
        return False, "Carbon connectivity graph is branched (expected 2 endpoints)"
    for node, nbrs in cgraph.items():
        if len(nbrs) > 2:
            return False, f"Carbon atom {node} has degree {len(nbrs)} (>2), indicating branching"
    
    # (C) The carboxyl carbon must be one of the endpoints.
    if carboxyl_carbon_idx not in endpoints:
        return False, "Carboxyl carbon is not at an end of the carbon chain"
    
    # Determine the ordered chain. Start from the carboxyl carbon and traverse the unique path.
    chain_order = [carboxyl_carbon_idx]
    prev = None
    current = carboxyl_carbon_idx
    while True:
        nbrs = cgraph[current]
        # Exclude the atom we just came from.
        next_atoms = [nbr for nbr in nbrs if nbr != prev]
        if not next_atoms:
            # reached the other endpoint
            break
        # There must be exactly one next atom.
        next_atom = next_atoms[0]
        chain_order.append(next_atom)
        prev, current = current, next_atom

    # Verify that the chain_order includes all carbons.
    if set(chain_order) != set(carbon_ids):
        return False, "Extra carbon atoms detected outside the main chain (branching present)"
    
    # For clarity, mark positions:
    chain_length = len(chain_order)
    # Identify the terminal indices:
    terminal_idxs = {chain_order[0], chain_order[-1]}  # one is carboxyl, the other is the methyl (or substituted terminal)
    
    # (E) Check substituents on each carbon in the chain.
    # Get a mapping from atom index to Atom object.
    atom_by_idx = {atom.GetIdx(): atom for atom in mol.GetAtoms()}
    for pos, c_idx in enumerate(chain_order):
        atom = atom_by_idx[c_idx]
        # Identify if this is carboxyl carbon.
        is_carboxyl = (c_idx == carboxyl_carbon_idx)
        is_terminal = (c_idx in terminal_idxs)
        # For each neighbor in the full molecule (including non-carbons)
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            # If neighbor is a carbon in the main chain then skip.
            if nbr_idx in chain_order:
                continue
            # Get the bond between atom and nbr.
            bond = mol.GetBondBetweenAtoms(c_idx, nbr_idx)
            # Allowed cases:
            if is_carboxyl:
                # Only allow oxygens that are part of the carboxyl group.
                if nbr.GetAtomicNum() != 8 or (nbr_idx not in allowed_for_carboxyl):
                    return False, f"Carboxyl carbon (atom {c_idx}) has an extra substituent (atom {nbr_idx})"
                continue  # allowed
            # For terminal carbon (that is not the carboxyl):
            if is_terminal:
                # Allow one -OH substituent: an oxygen attached by a SINGLE bond.
                # Note: If there is more than one non-chain substituent, disallow.
                if nbr.GetAtomicNum() == 8:
                    # If the bond is SINGLE then it is likely an -OH.
                    if bond.GetBondType() == Chem.BondType.SINGLE:
                        # Allowed: but if more than one such oxygen is attached, flag later.
                        pass
                    else:
                        return False, f"Terminal carbon (atom {c_idx}) has an improper oxygen bond (atom {nbr_idx})"
                else:
                    return False, f"Terminal carbon (atom {c_idx}) has an unexpected substituent (atom {nbr_idx}, atomic num {nbr.GetAtomicNum()})"
            else:
                # Internal carbon (not carboxyl and not terminal)
                # Allowed only if the substituent is an oxygen double-bonded (a keto group).
                if nbr.GetAtomicNum() == 8 and bond.GetBondType() == Chem.BondType.DOUBLE:
                    # Allowed keto substituent.
                    continue
                else:
                    return False, f"Internal chain carbon (atom {c_idx}) has a disallowed substituent (atom {nbr_idx})"
        # For terminal (non-carboxyl) carbon, ensure that if an oxygen is attached it is single and not multiple.
        if is_terminal and not is_carboxyl:
            oxygens = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(c_idx, nbr.GetIdx()).GetBondType() == Chem.BondType.SINGLE]
            if len(oxygens) > 1:
                return False, f"Terminal carbon (atom {c_idx}) has multiple oxygen substituents"
    
    return True, "Molecule is a straight‐chain saturated fatty acid with no side‐chain substituents"

# Example usage (for testing purposes)
if __name__ == "__main__":
    # We include a few examples that should be accepted as true positives.
    test_cases = [
        ("CC(=O)CCC(O)=O", "4-oxopentanoic acid"),
        ("CCCCCCCCCCCCCCCCCCCCCCCCCC(O)=O", "heptacosanoic acid"),
        ("CCCCCCCCCCCCCCCCCCCCCCCCCCCC(O)=O", "nonacosanoic acid"),
        ("CCCCCCCCCCCCCCCCCCCCCC(O)=O", "tricosanoic acid"),
        ("C(CCCCCCCCCCCCCCCC)CCCCCCCCCCC(O)=O", "octacosanoic acid"),
        ("CCCCCCCCCCCCCCCC(O)=O", "heptadecanoic acid"),
        ("CCCCCCCCCCCCCCCCCCCCCCCCC(O)=O", "hexacosanoic acid"),
        ("CCCCCCCCCCCCCCCCCCCCCCCCCCCCC(O)=O", "triacontanoic acid"),
        ("CCCCCCCCCCCCCCCC(O)=O", "hexadecanoic acid"),
        ("CCCCCCC(O)=O", "heptanoic acid"),
        ("OCCCCCCCCCCCCCCC(O)=O", "15-hydroxypentadecanoic acid"),
        ("CCCC(O)=O", "butyric acid"),
        ("CCCCCCCCCCCCCCCCCCCC(O)=O", "henicosanoic acid"),
        ("CCCCCCCCCCCCCCCCCCCCCCCCC(O)=O", "docosanoic acid"),
        ("CCCCCCCC(O)=O", "octanoic acid"),
        ("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC(O)=O", "dotriacontanoic acid"),
        ("CCCCCCCCCCCCCC(O)=O", "tetradecanoic acid"),
        ("OCCCCCCCCCCCCCCCCCCCC(O)=O", "20-hydroxyicosanoic acid"),
        ("CCCCCCCCCCCCCCCCCCC(O)=O", "nonadecanoic acid"),
        ("CCCCCC(O)=O", "hexanoic acid"),
        ("CCCCCCCCCCCCCCCCCCCCCC(O)=O", "icosanoic acid"),
        ("CCCCCCCCCCCC(O)=O", "tridecanoic acid"),
        ("OCCCCCCCCCC(O)=O", "10-hydroxycapric acid"),
        ("C(CCCCCC)CCC(=O)O", "decanoic acid"),
        ("OCCCCCCCCCCCCC(O)=O", "13-hydroxytridecanoic acid"),
    ]
    for smi, name in test_cases:
        valid, reason = is_straight_chain_saturated_fatty_acid(smi)
        print(f"SMILES: {smi}\n NAME: {name}\n -> {valid}: {reason}\n")