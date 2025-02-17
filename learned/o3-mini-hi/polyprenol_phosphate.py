"""
Classifies: CHEBI:16460 polyprenol phosphate
"""
"""
Classifies: Polyprenol phosphate
Definition: A prenol phosphate resulting from the formal condensation of the terminal 
          allylic (or isoprenoid) hydroxy group of a polyprenol with 1 mol eq. of phosphoric acid.

This implementation uses revised heuristics:
  1. For each phosphorus atom, find any oxygen neighbor that links to a candidate carbon.
  2. Only consider candidate carbons that are sp3‐hybridized, non‐aromatic and not in a ring.
  3. First, check if the candidate is immediately “allylic” – one of its other neighbors is involved in a non‐aromatic C=C double bond.
  4. If not, perform a DFS along acyclic, non‐aromatic, sp3 carbon atoms (the putative isoprenoid chain).
     The chain must be long (>=10 carbons) and contain several (>=3) double bonds.
  5. Also perform global checks on sufficient total carbons, presence of at least one non‐aromatic C=C bond, and a minimum molecular weight.
If these conditions are met, the function returns True with an explanation.
Otherwise, it returns False with a reason.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.rdchem import BondType

def is_polyprenol_phosphate(smiles: str):
    """
    Determines whether a molecule is a polyprenol phosphate.
    
    A polyprenol phosphate is defined as a prenol phosphate arising from the condensation
    of the terminal allylic (or isoprenoid) hydroxy group of a polyprenol with phosphoric acid.
    
    The function uses two sets of criteria:
      - Local criteria: from any phosphorus atom, an oxygen must connect to a candidate (sp3, acyclic)
        carbon that is either directly allylic (has a neighbor with a non‐aromatic C=C bond) or 
        is part of a long chain (>=10 carbons with >=3 C=C bonds along an acyclic, sp3 path).
      - Global criteria: the molecule overall must contain a high enough number of carbons, at least one
        non‐aromatic C=C bond and a molecular weight typical for polyprenol phosphates.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): True with an explanation if the molecule meets the criteria,
                     False with a reason otherwise.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    found_linkage = False
    linkage_reason = ""
    
    # Helper: Check if an atom is directly allylic.
    def is_directly_allylic(carbon):
        # Look at neighbors of candidate carbon (except the one coming from P–O)
        for nbr in carbon.GetNeighbors():
            if nbr.GetAtomicNum() != 6:
                continue
            # In our view an allylic situation is if any bond from this neighbor (besides the bond linking to candidate)
            # is a double bond (non‐aromatic) to another carbon.
            for bond in nbr.GetBonds():
                # Skip if this bond is between nbr and candidate
                if bond.GetBeginAtomIdx() == carbon.GetIdx() or bond.GetEndAtomIdx() == carbon.GetIdx():
                    continue
                # Only consider bonds between carbons that are non‐aromatic.
                other = bond.GetOtherAtom(nbr)
                if other.GetAtomicNum() != 6:
                    continue
                if not (nbr.GetIsAromatic() or other.GetIsAromatic()):
                    if bond.GetBondType() == BondType.DOUBLE:
                        return True
        return False

    # Helper: DFS traversal from a candidate carbon through acyclic, sp3, non‐aromatic carbons.
    def dfs_chain(start):
        visited = set()
        stack = [(start, None)]  # store tuples of (atom, bond from its parent)
        count_carbons = 0
        count_double = 0
        while stack:
            current, from_bond = stack.pop()
            if current.GetIdx() in visited:
                continue
            visited.add(current.GetIdx())
            count_carbons += 1
            if from_bond is not None and from_bond.GetBondType() == BondType.DOUBLE:
                count_double += 1
            for nb in current.GetNeighbors():
                # Only traverse to carbon neighbors that:
                # (a) are sp3 (so typically part of an isoprenoid chain),
                # (b) not aromatic, and (c) not in a ring.
                if nb.GetAtomicNum() == 6 and (nb.GetHybridization() == Chem.rdchem.HybridizationType.SP3) and (not nb.GetIsAromatic()) and (not nb.IsInRing()):
                    if nb.GetIdx() not in visited:
                        bond = mol.GetBondBetweenAtoms(current.GetIdx(), nb.GetIdx())
                        if bond is not None:
                            stack.append((nb, bond))
        return count_carbons, count_double

    # Now search among phosphorus atoms.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 15:
            continue  # not phosphorus
        # For each phosphorus atom, look at its oxygen neighbors.
        for o_atom in atom.GetNeighbors():
            if o_atom.GetAtomicNum() != 8:
                continue  # not oxygen
            # For each neighbor of the oxygen, find a candidate carbon.
            for cand in o_atom.GetNeighbors():
                if cand.GetIdx() == atom.GetIdx():
                    continue
                if cand.GetAtomicNum() != 6:
                    continue  # must be carbon
                # Require candidate to be sp3, non‐aromatic, and not in a ring.
                if cand.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
                    continue
                if cand.GetIsAromatic() or cand.IsInRing():
                    continue

                # First, check immediate (direct) allylicity.
                if is_directly_allylic(cand):
                    found_linkage = True
                    linkage_reason = "Allylic phosphate linkage found (directly attached allylic carbon)"
                    break

                # Otherwise, search along the chain from this candidate.
                chain_carbons, chain_double = dfs_chain(cand)
                if chain_carbons >= 10 and chain_double >= 3:
                    found_linkage = True
                    linkage_reason = (f"Phosphate linked to a long isoprenoid chain "
                                      f"(chain carbons: {chain_carbons}, double bonds: {chain_double})")
                    break
            if found_linkage:
                break
        if found_linkage:
            break

    if not found_linkage:
        return False, ("No appropriate phosphate linkage found (phosphate not connected via an oxygen to a non‐aromatic, acyclic, sp³ candidate carbon "
                       "with direct allylic environment or a sufficiently long isoprenoid chain).")

    # Global checks – must have enough carbon atoms, some non‐aromatic C=C bonds, and typical molecular weight.
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if total_carbons < 10:
        return False, f"Too few carbons in the molecule ({total_carbons}) for a polyprenol chain."
    
    total_double_bonds = 0
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        if (bond.GetBondType() == BondType.DOUBLE and 
            a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6 and 
            not a1.GetIsAromatic() and not a2.GetIsAromatic()):
            total_double_bonds += 1
    if total_double_bonds < 1:
        return False, f"Not enough C=C bonds ({total_double_bonds}) to support a polyprenol structure."
    
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if mw < 300:
        return False, f"Molecular weight ({mw:.1f} Da) seems too low for a polyprenol phosphate."
    
    return True, f"Contains polyprenol chain with phosphate linkage; {linkage_reason}."

# Example usage (for testing):
if __name__ == '__main__':
    # Test with one true positive example:
    test_smiles = "CC(C)=CCC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C/CC\\C(C)=C/CC\\C(C)=C/CC\\C(C)=C/CC\\C(C)=C/CC\\C(C)=C/CC\\C(C)=C/COP(O)(O)=O"
    result, reason = is_polyprenol_phosphate(test_smiles)
    print("Result:", result)
    print("Reason:", reason)