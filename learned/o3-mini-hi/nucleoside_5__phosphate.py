"""
Classifies: CHEBI:16701 nucleoside 5'-phosphate
"""
"""
Classifies: Nucleoside 5'-phosphate
Definition: A ribosyl or deoxyribosyl derivative of a pyrimidine or purine base 
in which C-5 of the ribose ring is mono-, di-, tri- or tetra-phosphorylated.
This heuristic algorithm searches for a five-membered sugar ring (four carbons and one oxygen),
a phosphate group attached to a CH2 substituent (as expected at C-5),
and an attached nucleobase (an aromatic nitrogen-containing ring).
"""

from rdkit import Chem

def is_nucleoside_5__phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside 5'-phosphate based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule matches the criteria for a nucleoside 5'-phosphate, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # We look for a five-membered ring that might represent the sugar.
    # A ribose (or deoxyribose) ring has 5 atoms with 4 carbons and one oxygen.
    ring_info = mol.GetRingInfo()
    sugar_ring_found = False
    sugar_ring_atoms = None
    for ring in ring_info.AtomRings():
        if len(ring) == 5:
            num_oxygens = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
            num_carbons = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
            if num_oxygens == 1 and num_carbons == 4:
                sugar_ring_found = True
                sugar_ring_atoms = set(ring)
                break
    if not sugar_ring_found:
        return False, "No ribose/deoxyribose (five-membered, 4C+1O) sugar ring found"
    
    # Next, check for an exocyclic CH2 group attached to the sugar ring that bears a phosphate.
    # We look for a carbon (with at least two hydrogens) bonded to a sugar ring atom
    # that, in turn, is connected (via an oxygen) to a phosphorus atom.
    phosphate_on_sugar = False
    for ring_idx in sugar_ring_atoms:
        sugar_atom = mol.GetAtomWithIdx(ring_idx)
        # Look at neighbors that are not part of the sugar ring
        for nbr in sugar_atom.GetNeighbors():
            if nbr.GetIdx() in sugar_ring_atoms:
                continue
            # Check if neighbor is a carbon likely being CH2
            if nbr.GetAtomicNum() == 6 and nbr.GetTotalNumHs() >= 2:
                # Now, check if this carbon has an oxygen neighbor that in turn is bonded to phosphorus.
                for sub_nbr in nbr.GetNeighbors():
                    if sub_nbr.GetAtomicNum() == 8:
                        # Check neighbors of the oxygen
                        for o_nbr in sub_nbr.GetNeighbors():
                            if o_nbr.GetAtomicNum() == 15:  # phosphorus atomic number
                                phosphate_on_sugar = True
                                break
                        if phosphate_on_sugar:
                            break
                if phosphate_on_sugar:
                    break
        if phosphate_on_sugar:
            break
    if not phosphate_on_sugar:
        return False, "No phosphate group found attached to a CH2 substituent of the sugar ring (expected at C-5)"
    
    # Lastly, search for an attached nucleobase.
    # We assume that a nucleobase is an aromatic heterocycle containing at least one nitrogen,
    # and is attached to one of the sugar ring atoms.
    base_found = False
    for ring_idx in sugar_ring_atoms:
        sugar_atom = mol.GetAtomWithIdx(ring_idx)
        for nbr in sugar_atom.GetNeighbors():
            if nbr.GetIdx() in sugar_ring_atoms:
                continue
            # Check if neighbor is aromatic and has at least one nitrogen in its immediate environment.
            if nbr.GetIsAromatic():
                # Check if this atom or its neighbors include nitrogen 
                if nbr.GetAtomicNum() == 7:
                    base_found = True
                    break
                else:
                    for nbr2 in nbr.GetNeighbors():
                        if nbr2.GetAtomicNum() == 7 and nbr2.GetIsAromatic():
                            base_found = True
                            break
            if base_found:
                break
        if base_found:
            break
    if not base_found:
        return False, "No nucleobase (aromatic N-containing group) found attached to the sugar ring"
    
    # If all criteria are met, return positive classification.
    return True, "Molecule contains a sugar ring with phosphate at C-5 and an attached nucleobase"

# (Optional) For testing, one might run:
# test_smiles = "O[C@@H]1[C@@H](COP(O)(O)=O)O[C@H]([C@@H]1O)n1ccc(=O)[nH]c1=O"  # uridine 5'-monophosphate
# result, reason = is_nucleoside_5__phosphate(test_smiles)
# print(result, reason)