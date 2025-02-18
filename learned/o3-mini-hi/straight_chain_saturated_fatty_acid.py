"""
Classifies: CHEBI:39418 straight-chain saturated fatty acid
"""
#!/usr/bin/env python
"""
Classifies: Straight-chain saturated fatty acid
Definition: Any saturated fatty acid lacking a side-chain.
Improved rules here:
  (1) Only single bonds between carbons (saturated).
  (2) Exactly one carboxyl group (COOH or COO–) is present.
  (3) When considering only carbon atoms connected by single bonds, they must
      form a single linear (unbranched) chain.
  (4) The carboxyl carbon must be at one end of that chain.
  (5) Off-chain substituents on carbons are allowed only as follows:
         - Carboxyl carbon: only the oxygens of the carboxyl group.
         - The other terminal carbon (methyl end): at most one oxygen (–OH) is allowed.
         - Internal carbons: must not have any non‐hydrogen substituents.
Any violation leads to rejection.
"""
from rdkit import Chem

def is_straight_chain_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a straight‐chain saturated fatty acid.
    Args:
         smiles (str): SMILES string of the molecule.
    Returns:
         (bool, str): Tuple where the boolean indicates whether the molecule meets
                      the criteria and the string explains the reasoning.
    """
    # Parse SMILES and add explicit hydrogens (to see substituents clearly).
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)
    
    # (1) Ensure all carbon–carbon bonds are single bonds.
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
            if bond.GetBondType() != Chem.BondType.SINGLE:
                return False, "Found one or more unsaturated C–C bonds."
    
    # (2) Find exactly one carboxyl group.
    # Use a SMARTS pattern that covers COOH and COO– groups.
    carboxyl_smarts = "C(=O)[O;H1,-1]"
    carboxyl_pattern = Chem.MolFromSmarts(carboxyl_smarts)
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    if len(carboxyl_matches) != 1:
        return False, f"Found {len(carboxyl_matches)} carboxyl groups; expected exactly 1."
    # In the match, assume the first atom is the carboxyl carbon and the subsequent atoms are the oxygens.
    carboxyl_match = carboxyl_matches[0]
    carboxyl_carbon = carboxyl_match[0]
    allowed_carboxyl_oxygens = set(carboxyl_match[1:])
    
    # (3) Build the connectivity graph among carbon atoms (only using single bonds).
    carbon_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    carbon_graph = { idx: [] for idx in carbon_indices }
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.SINGLE:
            i = bond.GetBeginAtomIdx()
            j = bond.GetEndAtomIdx()
            if i in carbon_graph and j in carbon_graph:
                carbon_graph[i].append(j)
                carbon_graph[j].append(i)
    # Ensure the carboxyl carbon is a carbon.
    if carboxyl_carbon not in carbon_graph:
        return False, "Carboxyl carbon not found among carbon atoms."
    
    # Use a DFS from the carboxyl carbon to see if all carbons belong to a single connected component.
    visited = set()
    stack = [carboxyl_carbon]
    while stack:
        node = stack.pop()
        if node not in visited:
            visited.add(node)
            for nbr in carbon_graph[node]:
                if nbr not in visited:
                    stack.append(nbr)
    if set(carbon_indices) != visited:
        return False, "Extra carbon atoms detected outside the main chain (branching present)."
    
    # (4) Check that the carbon connectivity graph forms a linear chain (path).
    # In a valid chain, exactly two carbons have degree 1 (the endpoints) and all others degree 2.
    endpoints = [node for node, nbrs in carbon_graph.items() if len(nbrs)==1]
    if len(endpoints) != 2:
        return False, "Carbon connectivity graph is branched (expected exactly 2 endpoints)."
    # The carboxyl carbon must be one of the endpoints.
    if carboxyl_carbon not in endpoints:
        return False, "Carboxyl carbon is not at one end of the carbon chain."
    
    # (5) Determine the linear chain order by traversing from the carboxyl carbon.
    chain_order = []
    current = carboxyl_carbon
    prev = None
    while True:
        chain_order.append(current)
        nbrs = [nbr for nbr in carbon_graph[current] if nbr != prev]
        if not nbrs:
            break
        if len(nbrs) != 1:
            return False, "Ambiguous carbon chain structure detected."
        prev, current = current, nbrs[0]
    
    # Make a lookup for atoms by index.
    atom_by_idx = {atom.GetIdx(): atom for atom in mol.GetAtoms()}
    # Identify the non-acid endpoint.
    other_endpoint = [pt for pt in endpoints if pt != carboxyl_carbon][0]
    
    # (6) For each carbon in the chain, examine all substituents (neighbors not in the chain).
    # Allowed:
    #  - At the carboxyl carbon: only oxygens that are part of the carboxyl group.
    #  - At the terminal non-carboxyl carbon (methyl end): allow at most one extra oxygen, and it must be in an –OH group.
    #  - For all internal carbons: no substituents (except hydrogens).
    for c_idx in chain_order:
        atom = atom_by_idx[c_idx]
        # Gather all substituents on this carbon that are NOT in the main chain.
        non_chain_subs = []
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() in chain_order:
                continue
            # Allow hydrogens (including deuterium which is atomic number 1).
            if nbr.GetAtomicNum() == 1:
                continue
            non_chain_subs.append(nbr)
        # Now check rules for each carbon:
        if c_idx == carboxyl_carbon:
            # Only allowed substituents are the carboxyl oxygens.
            for sub in non_chain_subs:
                if sub.GetIdx() not in allowed_carboxyl_oxygens:
                    return False, f"Carboxyl carbon (atom {c_idx}) has an extra substituent (atom {sub.GetIdx()})."
        elif c_idx == other_endpoint:
            # Terminal (non-carboxyl) carbon – allow at most one oxygen; if present it should be attached via a single bond.
            if len(non_chain_subs) > 1:
                return False, f"Terminal carbon (atom {c_idx}) has multiple substituents."
            if len(non_chain_subs) == 1:
                sub = non_chain_subs[0]
                bond = mol.GetBondBetweenAtoms(c_idx, sub.GetIdx())
                if sub.GetAtomicNum() != 8 or bond.GetBondType() != Chem.BondType.SINGLE:
                    return False, f"Terminal carbon (atom {c_idx}) has an unexpected substituent (atom {sub.GetIdx()}, atomic num {sub.GetAtomicNum()})."
        else:
            # Internal carbon: no substituents (other than hydrogens) permitted.
            if non_chain_subs:
                return False, f"Internal chain carbon (atom {c_idx}) has a disallowed substituent (atom {non_chain_subs[0].GetIdx()})."
    
    return True, "Molecule is a straight‐chain saturated fatty acid with no side‐chain substituents"

# Example usage:
if __name__ == "__main__":
    test_cases = [
        ("CC(=O)CCC(O)=O", "4-oxopentanoic acid"),
        ("CCCCCCCCCCCCCCCCCCCCCCCCCC(O)=O", "heptacosanoic acid"),
        ("CCCCCCCCCCCCCCCCCCCCCCCCCCCC(O)=O", "nonacosanoic acid"),
        ("CCCCCCCCCCCCCCCCCCCCCC(O)=O", "tricosanoic acid"),
        ("C(C(C(C(C(C(C(C(C(C(C(C(C(C(=O)O)([2H])[2H])([2H])[2H])([2H])[2H])([2H])[2H])([2H])[2H])([2H])[2H])([2H])[2H])([2H])[2H])([2H])[2H])([2H])[2H])([2H])([2H])[2H]", "tetradecanoic-d27 acid"),
        ("C(CCCCCCCCCCCCCCCC)CCCCCCCCCCC(O)=O", "octacosanoic acid"),
        ("CCCCCCCCCCCCCCCC(O)=O", "heptadecanoic acid"),
        ("CCCCCCCCCCCCCCCCCCCCCCCCC(O)=O", "hexacosanoic acid"),
        ("CCCCCCCCCCCCCCCCCCCCCCCCCCCCC(O)=O", "triacontanoic acid"),
        ("[2H]C([2H])([2H])C([2H])([2H])C([2H])([2H])C([2H])([2H])C([2H])([2H])C([2H])([2H])C([2H])([2H])C([2H])([2H])C([2H])([2H])C([2H])([2H])C([2H])([2H])C([2H])([2H])C(O)=O", "palmitic acid-d31"),
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
        # The following examples (with extra functional groups) should be rejected:
        ("OCCCCCC([O-])=O", "6-hydroxyhexanoate"),
        ("CCCCCCCC([O-])=O", "nonanoate"),
    ]
    
    for smi, name in test_cases:
        valid, reason = is_straight_chain_saturated_fatty_acid(smi)
        print(f"SMILES: {smi}\n NAME: {name}\n  -> {valid}: {reason}\n")