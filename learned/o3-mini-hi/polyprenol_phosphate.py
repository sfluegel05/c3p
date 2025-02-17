"""
Classifies: CHEBI:16460 polyprenol phosphate
"""
"""
Classifies: Polyprenol phosphate
Definition: A prenol phosphate resulting from the formal condensation of the terminal
          allylic (or isoprenoid) hydroxy group of a polyprenol with 1 mol eq. of phosphoric acid.

This implementation uses the following heuristics:
  1. Find at least one phosphorus atom whose oxygen neighbor is connected to a carbon.
  2. Require that the attached carbon is not aromatic and not in a ring.
  3. Check if the candidate carbon is “allylic” – i.e. one of its other neighboring carbons is attached to a C=C double bond.
  4. If not directly allylic, follow (via DFS) the candidate carbon along an acyclic, non‐aromatic carbon chain.
     The chain must be long enough (>=10 carbon atoms) and contain at least 3 C=C bonds.
  5. Also, require that the molecule globally has a sufficient number of carbon atoms,
     some C=C bonds (and a molecular weight typical of polyprenol phosphates).
If these conditions are met, the function returns True with an explanation.
Otherwise it returns False with a reason.
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
        str: Explanation of the decision.
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    found_linkage = False
    linkage_reason = ""
    
    # We look through each phosphorus atom.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 15:
            continue  # not phosphorus
        # For each phosphorus, check its oxygen neighbors
        for o_atom in atom.GetNeighbors():
            if o_atom.GetAtomicNum() != 8:
                continue  # not oxygen
            # Now check each neighbor of the oxygen (candidate carbon)
            for cand in o_atom.GetNeighbors():
                if cand.GetIdx() == atom.GetIdx():
                    continue  # skip back‐bond to phosphorus
                if cand.GetAtomicNum() != 6:
                    continue  # must be carbon
                # Only consider candidates that are non‐aromatic and not in a ring.
                if cand.GetIsAromatic() or cand.IsInRing():
                    continue
                
                # First heuristic: check if candidate carbon is directly allylic.
                immediate_allylic = False
                for nbr in cand.GetNeighbors():
                    if nbr.GetIdx() == o_atom.GetIdx():
                        continue  # skip the oxygen we came from
                    if nbr.GetAtomicNum() != 6:
                        continue  # only consider carbons
                    # Check all bonds of this neighbor for a C=C double bond (and that the bond is not aromatic)
                    for bond in nbr.GetBonds():
                        # Get the other atom in the bond:
                        begin_idx = bond.GetBeginAtomIdx()
                        end_idx = bond.GetEndAtomIdx()
                        other_idx = begin_idx if end_idx == nbr.GetIdx() else end_idx
                        if other_idx == cand.GetIdx():
                            continue
                        other_atom = mol.GetAtomWithIdx(other_idx)
                        # Only consider non‐aromatic C=C bonds
                        if other_atom.GetAtomicNum() == 6 and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                            immediate_allylic = True
                            break
                    if immediate_allylic:
                        break
                if immediate_allylic:
                    found_linkage = True
                    linkage_reason = "Allylic phosphate linkage found (directly attached allylic carbon)"
                    break  # we found a valid linkage
                
                # Second heuristic: if not immediately allylic check if the candidate is connected to a long isoprenoid chain.
                # Traverse along carbon atoms (ignoring rings and aromatic atoms) starting from candidate.
                visited = set()
                stack = [(cand, None)]  # (current atom, parent bond type info)
                chain_carbons = 0
                chain_double_bonds = 0
                # We perform DFS over acyclic, non‐aromatic carbons.
                while stack:
                    current, via_bond = stack.pop()
                    if current.GetIdx() in visited:
                        continue
                    visited.add(current.GetIdx())
                    # Count current carbon.
                    chain_carbons += 1
                    # If we reached this atom via a bond, count if that bond is a double bond.
                    if via_bond is not None and via_bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        chain_double_bonds += 1
                    # Expand along neighbors that are carbons, non‐aromatic and not in rings.
                    for nb in current.GetNeighbors():
                        if nb.GetAtomicNum() == 6 and (not nb.GetIsAromatic()) and (not nb.IsInRing()):
                            # Avoid going backwards to the oxygen that attached to phosphorus:
                            if nb.GetIdx() not in visited:
                                bond = mol.GetBondBetweenAtoms(current.GetIdx(), nb.GetIdx())
                                if bond is not None:
                                    stack.append((nb, bond))
                if chain_carbons >= 10 and chain_double_bonds >= 3:
                    found_linkage = True
                    linkage_reason = (f"Phosphate linked to a long isoprenoid chain "
                                      f"(chain carbons: {chain_carbons}, double bonds: {chain_double_bonds})")
                    break
            if found_linkage:
                break
        if found_linkage:
            break
    
    if not found_linkage:
        return False, "No appropriate phosphate linkage found (phosphate not connected via an oxygen to a non‐aromatic, acyclic candidate carbon with typical allylic/isoprenoid features)"
    
    # Global checks:
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if total_carbons < 10:
        return False, f"Too few carbons in the molecule ({total_carbons}) for a polyprenol chain."
    
    total_double_bonds = 0
    for bond in mol.GetBonds():
        # Only count non‐aromatic C=C bonds between carbons.
        if (bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and 
            bond.GetBeginAtom().GetAtomicNum() == 6 and 
            bond.GetEndAtom().GetAtomicNum() == 6 and
            not bond.GetBeginAtom().GetIsAromatic() and
            not bond.GetEndAtom().GetIsAromatic()):
            total_double_bonds += 1
    if total_double_bonds < 1:
        return False, f"Not enough C=C bonds ({total_double_bonds}) to support a polyprenol structure."
    
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if mw < 300:
        return False, f"Molecular weight ({mw:.1f} Da) seems too low for a polyprenol phosphate."
    
    return True, f"Contains polyprenol chain with phosphate linkage; {linkage_reason}."

# Example usage (for testing):
if __name__ == '__main__':
    # Test with one of the true positive examples:
    test_smiles = "CC(C)=CCC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C/CC\\C(C)=C/CC\\C(C)=C/CC\\C(C)=C/CC\\C(C)=C/CC\\C(C)=C/CC\\C(C)=C/COP(O)(O)=O"
    result, reason = is_polyprenol_phosphate(test_smiles)
    print("Result:", result)
    print("Reason:", reason)