"""
Classifies: CHEBI:16460 polyprenol phosphate
"""
"""
Classifies: Polyprenol phosphate
Definition: A prenol phosphate resulting from the formal condensation of 
the terminal allylic hydroxy group of a polyprenol with 1 mol eq. of phosphoric acid.

This implementation uses the following heuristics:
  1. It searches for a phosphate group (P) that is connected via an oxygen to a carbon.
  2. For each such candidate carbon it first checks for an allylic linkage 
     (i.e. at least one neighboring carbon is involved in a C=C double bond).
  3. If the allylic check fails then it “follows” the candidate carbon along a chain
     (ignoring ring atoms) to see if the chain is long (>=10 carbon atoms) and contains 
     at least one C=C bond. This allows for cases where the terminal is saturated or masked.
  4. It also requires (globally) that the molecule has a minimal number of carbon atoms 
     and at least one C=C bond, and that its molecular weight is within a typical range.
If these conditions are met, the function returns True with an explanation. Otherwise False.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polyprenol_phosphate(smiles: str):
    """
    Determines if a molecule belongs to the polyprenol phosphate class.
    A polyprenol phosphate is defined as a prenol phosphate arising from the condensation
    of the terminal allylic (or isoprenoid) hydroxy group of a polyprenol with phosphoric acid.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is a polyprenol phosphate, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Initialize flag and reason for finding an appropriate phosphate linkage.
    found_linkage = False
    linkage_reason = ""
    
    # Look for a phosphorus atom.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 15:  # only phosphorus atoms
            continue
        # For each phosphorus atom, examine its oxygen neighbors.
        for o_atom in atom.GetNeighbors():
            if o_atom.GetAtomicNum() != 8:
                continue
            # Now check for a candidate carbon attached to the oxygen.
            for cand in o_atom.GetNeighbors():
                if cand.GetIdx() == atom.GetIdx():
                    continue  # skip phosphorus back‐link
                if cand.GetAtomicNum() != 6:
                    continue  # candidate must be carbon
                # First try: is candidate carbon directly “allylic”?
                allylic = False
                for nbr in cand.GetNeighbors():
                    if nbr.GetIdx() == o_atom.GetIdx():
                        continue  # skip the O we came from
                    if nbr.GetAtomicNum() != 6:
                        continue
                    # Check all bonds of this neighbor for a carbon–carbon double bond.
                    for bond in nbr.GetBonds():
                        # Get the other atom in the bond
                        begin_idx = bond.GetBeginAtomIdx()
                        end_idx = bond.GetEndAtomIdx()
                        other_idx = begin_idx if end_idx == nbr.GetIdx() else end_idx
                        if other_idx == cand.GetIdx():
                            continue
                        other_atom = mol.GetAtomWithIdx(other_idx)
                        if other_atom.GetAtomicNum() != 6:
                            continue
                        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                            allylic = True
                            break
                    if allylic:
                        break
                if allylic:
                    found_linkage = True
                    linkage_reason = "Allylic phosphate linkage found"
                    break  # break out: we found a valid connection.
                else:
                    # Second try: if the candidate is not obviously allylic, check the chain.
                    # Only consider candidate carbons that are not in rings.
                    if not cand.IsInRing():
                        visited = set()
                        stack = [cand]
                        carbon_count = 0
                        double_bonds_in_chain = 0
                        # Traverse (DFS) along carbon atoms, ignoring rings.
                        while stack:
                            current = stack.pop()
                            if current.GetIdx() in visited:
                                continue
                            visited.add(current.GetIdx())
                            if current.GetAtomicNum() == 6:
                                carbon_count += 1
                            for neigh in current.GetNeighbors():
                                if neigh.GetAtomicNum() == 6 and not neigh.IsInRing():
                                    bond = mol.GetBondBetweenAtoms(current.GetIdx(), neigh.GetIdx())
                                    if bond is not None:
                                        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                                            double_bonds_in_chain += 1
                                        if neigh.GetIdx() not in visited:
                                            stack.append(neigh)
                        # Consider long isoprenoid chains even if the terminal carbon is not strictly allylic.
                        if carbon_count >= 10 and double_bonds_in_chain >= 1:
                            found_linkage = True
                            linkage_reason = (f"Phosphate linked to a long isoprenoid chain "
                                               f"(chain carbons: {carbon_count}, "
                                               f"chain double bonds: {double_bonds_in_chain})")
                            break
            if found_linkage:
                break
        if found_linkage:
            break
    
    if not found_linkage:
        return False, "No appropriate phosphate linkage found (phosphate not connected via an oxygen to a candidate isoprenoid carbon)"
    
    # Global checks: many polyprenols have multiple C=C bonds and a long carbon chain.
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if total_carbons < 10:
        return False, f"Too few carbon atoms ({total_carbons}) for a polyprenol chain."
    
    total_double_bonds = 0
    for bond in mol.GetBonds():
        if (bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and 
            bond.GetBeginAtom().GetAtomicNum() == 6 and 
            bond.GetEndAtom().GetAtomicNum() == 6):
            total_double_bonds += 1
    if total_double_bonds < 1:
        return False, f"Not enough C=C bonds ({total_double_bonds}) to support a polyprenol structure."
    
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if mw < 300:
        return False, f"Molecular weight ({mw:.1f} Da) seems too low for a polyprenol phosphate."
    
    return True, f"Contains polyprenol chain with phosphate linkage; {linkage_reason}."

# Example usage (for testing):
if __name__ == '__main__':
    # Test with one of the provided true-positive examples:
    test_smiles = "CC(C)=CCC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C/CC\\C(C)=C/CC\\C(C)=C/CC\\C(C)=C/CC\\C(C)=C/CC\\C(C)=C/CC\\C(C)=C/COP(O)(O)=O"
    result, reason = is_polyprenol_phosphate(test_smiles)
    print("Result:", result)
    print("Reason:", reason)