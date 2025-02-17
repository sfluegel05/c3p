"""
Classifies: CHEBI:59238 cyclic fatty acid
"""
"""
Classifies: Cyclic Fatty Acid
Definition: Any fatty acid containing anywhere in its structure a ring of atoms.
A cyclic fatty acid here is defined (for our purposes) as a molecule that:
  (1) contains at least one ring,
  (2) does not contain amide bonds (to avoid peptides),
  (3) contains a terminal fatty acyl group – that is, a carbonyl functional group (acid or ester)
      with an oxygen (–OH or –OR) that is not attached to any nitrogen and is connected to exactly one
      acyclic, non‐aromatic carbon, and
  (4) from that attachment point there is a contiguous unbranched acyclic (non‐aromatic) carbon “tail” (chain)
      of at least 5 atoms.
      
If no valid fatty acyl group is detected under these rules, the function returns False.
"""

from rdkit import Chem

def is_cyclic_fatty_acid(smiles: str):
    """
    Determines whether the given SMILES string represents a cyclic fatty acid.
    
    Steps:
      1. Parse the SMILES and ensure the molecule exists.
      2. Check that at least one ring is present.
      3. Reject molecules containing amide bonds (SMARTS: [NX3][CX3](=O)).
      4. Look for candidate fatty acyl groups:
         - Find a carbon (atomic number 6) that is double bonded to an oxygen.
         - Ensure that this carbon is not attached to any nitrogen.
         - From its neighbors, select those that are not oxygen.
         - Require that exactly one such non‐oxygen neighbor exists.
         - This neighbor becomes the start of the fatty acyl "tail."
      5. Confirm that the attachment (tail starting carbon) is an acyclic, non‐aromatic carbon.
      6. Create a subgraph of acyclic, non‐aromatic carbons (by iterating over atoms and selecting carbons
         that are not in any ring and are not aromatic).
      7. Within this subgraph, verify that the start atom is terminal (degree 1) and “walk” linearly down the chain.
         Count the number of carbons in this contiguous chain. A valid chain has at least 5 carbons.
         
    Args:
        smiles (str): SMILES string representing the molecule.
        
    Returns:
        bool: True if the molecule meets the criteria for a cyclic fatty acid.
        str: Explanation of the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Check: molecule must have at least one ring.
    if not mol.GetRingInfo().AtomRings():
        return False, "No ring found in the structure"
    
    # 2. Reject molecule if any amide bonds present.
    # Pattern for amide: nitrogen-single-bonded to a carbonyl.
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=O)")
    if mol.HasSubstructMatch(amide_pattern):
        return False, "Molecule contains amide bonds (likely a peptide)"
    
    valid_chain_lengths = []
    
    # 3. Iterate over all atoms to find candidate carbonyl carbons.
    for atom in mol.GetAtoms():
        # Look only at carbon atoms.
        if atom.GetAtomicNum() != 6:
            continue

        # Check for a double bond to oxygen (i.e. a carbonyl group).
        # We loop over bonds and look for a double bond with an oxygen.
        has_carbonyl = False
        oxygen_partner = None
        for bond in atom.GetBonds():
            # Check if bond is double and one end is oxygen.
            if bond.GetBondTypeAsDouble() == 2:
                nbr = bond.GetOtherAtom(atom)
                if nbr.GetAtomicNum() == 8:
                    has_carbonyl = True
                    oxygen_partner = nbr
                    break
        if not has_carbonyl:
            continue
        
        # Exclude candidate if the carbonyl carbon is attached to any nitrogen (to avoid amide scenario).
        if any(nbr.GetAtomicNum() == 7 for nbr in atom.GetNeighbors()):
            continue
        
        # Identify non-oxygen neighbors (should be exactly one, the acyl attachment point).
        non_oxygen_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() != 8]
        if len(non_oxygen_neighbors) != 1:
            continue
        
        start_atom = non_oxygen_neighbors[0]
        # The attached atom must be carbon and should be acyclic and non‐aromatic.
        if start_atom.GetAtomicNum() != 6:
            continue
        if start_atom.IsInRing() or start_atom.GetIsAromatic():
            continue
        
        # 4. Build a set of candidate atoms representing acyclic (non‐ring, non‐aromatic) carbons.
        acyclic_atoms = {at.GetIdx() for at in mol.GetAtoms() 
                         if at.GetAtomicNum() == 6 and (not at.IsInRing()) and (not at.GetIsAromatic())}
        if start_atom.GetIdx() not in acyclic_atoms:
            continue  # Should not happen but safeguard.
        
        # Build an adjacency dictionary for the acyclic subgraph.
        neighbors_dict = {idx: [] for idx in acyclic_atoms}
        for idx in acyclic_atoms:
            a = mol.GetAtomWithIdx(idx)
            for nbr in a.GetNeighbors():
                if nbr.GetIdx() in acyclic_atoms:
                    neighbors_dict[idx].append(nbr.GetIdx())
        
        # For an unbranched (terminal) chain, the starting atom must have degree exactly 1.
        if len(neighbors_dict[start_atom.GetIdx()]) != 1:
            continue
        
        # 5. Walk linearly down the chain.
        chain_length = 1  # starting atom counts as one
        current = start_atom.GetIdx()
        previous = None
        while True:
            # Get acyclic neighbors excluding where we came from.
            next_atoms = [n for n in neighbors_dict[current] if n != previous]
            if len(next_atoms) == 1:
                chain_length += 1
                previous = current
                current = next_atoms[0]
            else:
                break
        
        valid_chain_lengths.append(chain_length)
    
    if not valid_chain_lengths:
        return False, "No valid terminal fatty acyl group found (acid/ester with a single attachment to an acyclic carbon)"
    
    max_chain = max(valid_chain_lengths)
    if max_chain < 5:
        return False, f"Found a terminal fatty acyl chain, but chain length is insufficient ({max_chain} carbons; need at least 5)"
    
    return True, f"Contains a terminal fatty acyl group with a contiguous acyclic carbon chain of length {max_chain} and a ring in the structure"

# For testing:
if __name__ == "__main__":
    test_examples = [
        # Examples (with provided example SMILES strings):
        ("O1C(CCCCCCCCCCC(OC[C@@H](OC(=O)CCCCCCCCCCC=2OC(CCCCC)=CC2C)CO)=O)=C(C(=C1CCCCC)C", "DG(11D5/11M5/0:0)"),
        ("OC(=O)CCCCCCNC1c2ccccc2CCc2ccccc12", "amineptine"),
        ("C(CCC/C=C\\C[C@@H]1/C(/O1)=C/C=C\\C/C=C\\CCCCC)(=O)O", "(5Z,8R,9Z,11Z,14Z)-8,9-epoxyicosatetraenoic acid"),
        ("S(C1CCC(OC(CCC(=O)C1)C)=O)CC(O)C(O)=O", "2-Hydroxy-3-(10-methyl-2,7-dioxooxecan-5-yl)sulfanylpropanoic acid"),
        ("O=C1OCC(=C1C)CCCC(=O)O", "Ficifuranone A"),
        ("O1C(CCCCC(O)=O)=C(C=C1CCCCC)", "3-Methyl-5-pentyl-2-furanpentanoic acid"),
        ("O1C(=CC=C1C(=O)C#C/C=C/C)/C=C\\C(O)=O", "Wyeronic acid"),
        ("OC(=O)CCCC[C@@H]1CCSS1", "(R)-lipoic acid")
    ]
    
    for smi, name in test_examples:
        result, reason = is_cyclic_fatty_acid(smi)
        print(f"Name: {name}\nSMILES: {smi}\nResult: {result}\nReason: {reason}\n{'-'*60}")