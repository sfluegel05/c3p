"""
Classifies: CHEBI:26244 prenols
"""
"""
Classifies: Prenols – any alcohol possessing the general formula 
    H-[CH2C(Me)=CHCH2]nOH
where the carbon skeleton is exclusively composed of one or more isoprene units.
This version uses global heuristics:
  - The molecule must be valid, acyclic, and contain only C and O.
  - The total number of carbons must be a multiple of 5.
  - It must have at least one free hydroxyl group (–OH not part of a carboxyl)
    that is attached to a terminal carbon (in the main carbon scaffold). For a diol,
    the two –OH groups must be at opposite termini.
  - The number of C=C bonds (between carbons) should equal the number of isoprene units,
    which is total carbons//5.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_prenols(smiles: str):
    """
    Determines if a molecule is a prenol (or prenol derivative) based on its SMILES.
    
    A prenol is defined as an alcohol whose carbon skeleton is exclusively
    composed of one or more isoprene units. In our heuristic, the molecule:
      - must be acyclic,
      - can contain only carbon and oxygen as heavy atoms,
      - must have a total number of carbons divisible by 5,
      - must have one (or two) free hydroxyl group(s) that are terminal (attached
        to a primary carbon of the longest carbon scaffold),
      - and must have a number of C=C bonds (between carbons) equal to (total_C/5).
    
    Args:
      smiles (str): SMILES string for the molecule.
      
    Returns:
      bool: True if classified as prenol, False otherwise.
      str: Explanation of the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Molecule must be acyclic.
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule is cyclic, not a prenol"

    # 2. Only allow carbon and oxygen as heavy atoms.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in (6, 8):
            return False, f"Molecule contains atom {atom.GetSymbol()} not in (C, O)"
    
    # 3. Count total number of carbon atoms.
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    carbon_count = len(carbon_atoms)
    if carbon_count < 5:
        return False, f"Too few carbons ({carbon_count}) to be a prenol"
    if carbon_count % 5 != 0:
        return False, f"Total carbon count ({carbon_count}) is not a multiple of 5; does not match isoprene repeats"
    expected_units = carbon_count // 5

    # 4. Identify free hydroxyl groups.
    # First, find all hydroxyl groups: pattern [OX2H]
    hydroxyl_pat = Chem.MolFromSmarts("[OX2H]")
    all_oh = mol.GetSubstructMatches(hydroxyl_pat)
    if not all_oh:
        return False, "No hydroxyl (–OH) group found"
    
    # Exclude hydroxyls that are part of carboxy groups: C(=O)[OX2H]
    carboxyl_pat = Chem.MolFromSmarts("C(=O)[OX2H]")
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pat)
    carboxyl_oxygens = set()
    for match in carboxyl_matches:
        # In the pattern, the oxygen in the hydroxyl is at position 1.
        carboxyl_oxygens.add(match[1])
    
    # Filter free hydroxyl groups (their oxygen atom index not in carboxyl_oxygens)
    free_oh = [match for match in all_oh if match[0] not in carboxyl_oxygens]
    if not free_oh:
        return False, "Only carboxylic acid hydroxyl(s) found; no free –OH present"
    
    # 5. Check that free –OH groups are terminal:
    # A hydroxyl is deemed “terminal” if its oxygen is attached to a primary carbon.
    # (i.e. the attached carbon should have exactly one heavy neighbor besides the –OH)
    terminal_oh_idxs = []
    for match in free_oh:
        o_atom = mol.GetAtomWithIdx(match[0])
        # Find heavy neighbors (ignore hydrogens)
        heavy_neighbors = [nbr for nbr in o_atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
        if not heavy_neighbors:
            continue  # should not happen for –OH though
        attached_c = heavy_neighbors[0]
        # Count heavy neighbors of the carbon (excluding the O we came from)
        c_heavy_neighbors = [nbr for nbr in attached_c.GetNeighbors() if nbr.GetAtomicNum() > 1 and nbr.GetIdx() != o_atom.GetIdx()]
        # For a primary carbon in a linear chain, expect 1 heavy neighbor.
        if len(c_heavy_neighbors) == 1:
            terminal_oh_idxs.append(o_atom.GetIdx())
    if len(terminal_oh_idxs) < 1:
        return False, "No terminal (primary) free –OH group found"
    if len(terminal_oh_idxs) > 2:
        return False, f"Too many terminal free –OH groups found ({len(terminal_oh_idxs)}); expected 1 or 2"
    
    # If two –OH groups are present, ensure they are at opposite ends.
    if len(terminal_oh_idxs) == 2:
        # Compute the attached carbon indices for each OH.
        attached_cs = []
        for o_idx in terminal_oh_idxs:
            o_atom = mol.GetAtomWithIdx(o_idx)
            for nbr in o_atom.GetNeighbors():
                if nbr.GetAtomicNum() == 6:
                    attached_cs.append(nbr.GetIdx())
                    break
        # Compute the shortest path between the two attached carbons.
        if len(attached_cs) == 2:
            path = Chem.GetShortestPath(mol, attached_cs[0], attached_cs[1])
            # In a linear prenol, the number of carbons in the longest chain should equal carbon_count.
            if len(path) != carbon_count:
                return False, "Terminal –OH groups are not at opposite ends of the carbon skeleton"
    
    # 6. Count carbon–carbon double bonds.
    double_bonds = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
                double_bonds += 1
    if double_bonds != expected_units:
        return False, (f"Mismatch in double bond count: found {double_bonds} C=C bond(s) "
                       f"but expected {expected_units} (from {carbon_count} carbons)")
    
    # If all tests pass, classify as prenol.
    # Count free hydroxyl groups (could be more than terminal in some borderline cases)
    return True, (f"Classified as prenol: contains {expected_units} isoprene unit(s) "
                  f"(by {carbon_count} carbons and {double_bonds} C=C bonds), "
                  f"{len(free_oh)} free –OH group(s) (with {len(terminal_oh_idxs)} terminal).")

# Example usage (for testing):
if __name__ == "__main__":
    # Example: prenol SMILES "CC(C)=CCO"
    test_smiles = "CC(C)=CCO"
    result, reason = is_prenols(test_smiles)
    print("Result:", result)
    print("Reason:", reason)