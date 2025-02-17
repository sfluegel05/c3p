"""
Classifies: CHEBI:39418 straight-chain saturated fatty acid
"""
#!/usr/bin/env python
"""
Classifies: Straight-chain saturated fatty acid
Definition: Any saturated fatty acid lacking a side-chain.
Rules implemented here:
  (1) Molecule must contain exactly one carboxyl group (COOH or COO–).
  (2) When considering only carbon atoms connected via single bonds,
      they must form precisely one linear (unbranched) chain.
  (3) The carboxyl carbon must be at one end of that chain.
  (4) For other carbons:
         - Terminal (non‐carboxyl) carbons: at most one extra -OH group (oxygen attached via a SINGLE bond)
         - Internal chain carbon: only a double-bonded oxygen (a keto oxygen) is allowed.
  (5) No unsaturation is allowed in the main chain (all C–C bonds must be single).
Any violation leads to a rejection.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_straight_chain_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a straight-chain saturated fatty acid.
    Args:
         smiles (str): SMILES string of the molecule.
    Returns:
         (bool, str): Tuple where the boolean indicates whether the molecule meets
                      the criteria and the string explains the reason.
    """
    # Parse SMILES and add explicit hydrogens (to catch substituents such as –OH)
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)
    
    # (1) Check that there are only single bonds between connected carbons.
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
            if bond.GetBondType() != Chem.BondType.SINGLE:
                return False, "Found an unsaturated carbon–carbon bond"
    
    # (2) Identify exactly one carboxyl group.
    # Using a SMARTS pattern for carboxyl groups (works for COOH and COO–)
    carboxyl_smarts = "C(=O)[O;H1,-1]"
    carboxyl_pattern = Chem.MolFromSmarts(carboxyl_smarts)
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    if len(carboxyl_matches) != 1:
        return False, f"Found {len(carboxyl_matches)} carboxyl groups; expected exactly 1"
    # In the matched tuple, assume the first atom is the carboxyl carbon.
    carboxyl_match = carboxyl_matches[0]
    carboxyl_carbon_idx = carboxyl_match[0]
    allowed_for_carboxyl = set(carboxyl_match[1:])  # Allowed oxygen atoms on the carboxyl carbon
    
    # (3) Build the connectivity graph for carbon atoms using only single bonds.
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    carbon_ids = [atom.GetIdx() for atom in carbon_atoms]
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
    
    # (4) Check that the carbon graph forms a linear chain:
    # Exactly two carbons should have degree 1 (endpoints) and all others degree 2.
    endpoints = [node for node, nbrs in cgraph.items() if len(nbrs)==1]
    if len(endpoints) != 2:
        return False, "Carbon connectivity graph is branched (expected 2 endpoints)"
    for node, nbrs in cgraph.items():
        if len(nbrs) > 2:
            return False, f"Carbon atom {node} has degree {len(nbrs)} (>2), indicating branching"
    
    # (5) The carboxyl carbon must be one of the endpoints.
    if carboxyl_carbon_idx not in endpoints:
        return False, "Carboxyl carbon is not at an end of the carbon chain"
    
    # (6) Determine the linear chain order by traversing from the carboxyl carbon.
    chain_order = [carboxyl_carbon_idx]
    prev = None
    current = carboxyl_carbon_idx
    while True:
        nbrs = cgraph[current]
        next_atoms = [nbr for nbr in nbrs if nbr != prev]
        if not next_atoms:
            break  # reached the far endpoint
        if len(next_atoms) != 1:
            return False, "Ambiguous path in carbon chain detected"
        next_atom = next_atoms[0]
        chain_order.append(next_atom)
        prev, current = current, next_atom
    if set(chain_order) != set(carbon_ids):
        return False, "Extra carbon atoms detected outside the main chain (branching present)"
    
    # Identify the two endpoints: carboxyl (must be one endpoint) and the other terminal.
    terminal_idxs = {chain_order[0], chain_order[-1]}
    
    # (7) Check substituents on each carbon in the chain.
    # Map atom index to atom for easy lookup.
    atom_by_idx = {atom.GetIdx(): atom for atom in mol.GetAtoms()}
    for pos, c_idx in enumerate(chain_order):
        atom = atom_by_idx[c_idx]
        is_carboxyl = (c_idx == carboxyl_carbon_idx)
        is_terminal = (c_idx in terminal_idxs)
        # Loop through neighbors (these include hydrogens and other substituents)
        # Ignore hydrogens (atomic num 1) as they are normal.
        # Also ignore neighbors which are part of the main carbon chain.
        extra_oxygens = []  # keep track of oxygen substituents attached to terminal (non-carboxyl)
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            # Ignore hydrogen substituents.
            if nbr.GetAtomicNum() == 1:
                continue
            # Skip if the neighbor is in the main chain.
            if nbr_idx in chain_order:
                continue
            # Get the bond between atom and nbr.
            bond = mol.GetBondBetweenAtoms(c_idx, nbr_idx)
            # (A) If this is the carboxyl carbon then only permit the oxygen atoms that are part of the carboxyl group.
            if is_carboxyl:
                if nbr.GetAtomicNum() != 8 or (nbr_idx not in allowed_for_carboxyl):
                    return False, f"Carboxyl carbon (atom {c_idx}) has an extra substituent (atom {nbr_idx})"
                continue
            # (B) For a terminal (non-carboxyl) carbon:
            if is_terminal:
                if nbr.GetAtomicNum() == 8:
                    # Allow an -OH if attached by a SINGLE bond.
                    if bond.GetBondType() == Chem.BondType.SINGLE:
                        extra_oxygens.append(nbr_idx)
                    else:
                        return False, f"Terminal carbon (atom {c_idx}) has an oxygen with improper bond (atom {nbr_idx})"
                else:
                    return False, f"Terminal carbon (atom {c_idx}) has an unexpected substituent (atom {nbr_idx}, atomic num {nbr.GetAtomicNum()})"
            else:
                # (C) For internal carbons: only allow a double-bonded oxygen (keto group).
                if nbr.GetAtomicNum() == 8 and bond.GetBondType() == Chem.BondType.DOUBLE:
                    continue
                else:
                    return False, f"Internal chain carbon (atom {c_idx}) has a disallowed substituent (atom {nbr_idx})"
        if is_terminal and (not is_carboxyl):
            if len(extra_oxygens) > 1:
                return False, f"Terminal carbon (atom {c_idx}) has multiple oxygen substituents"
    
    return True, "Molecule is a straight‐chain saturated fatty acid with no side‐chain substituents"

# Example usage for testing:
if __name__ == "__main__":
    test_cases = [
        ("CC(=O)CCC(O)=O", "4-oxopentanoic acid"),
        ("CCCCCCCCCCCCCCCCCCCCCCCCCC(O)=O", "heptacosanoic acid"),
        ("CCCCCCCCCCCCCCCCCCCCCCCCCCCC(O)=O", "nonacosanoic acid"),
        ("CCCCCCCCCCCCCCCCCCCCCC(O)=O", "tricosanoic acid"),
        ("C(C(CCCCCCCCCCCCCCCC)=O)CCCCCCCCCCC(O)=O", "octacosanoic acid"),
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
        # Some of the provided examples with deuterium or stereochemistry should also work.
    ]
    
    for smi, name in test_cases:
        valid, reason = is_straight_chain_saturated_fatty_acid(smi)
        print(f"SMILES: {smi}\n NAME: {name}\n  -> {valid}: {reason}\n")