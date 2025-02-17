"""
Classifies: CHEBI:16460 polyprenol phosphate
"""
"""
Classifies: Polyprenol phosphate
Definition: A prenol phosphate resulting from the formal condensation of 
the terminal allylic hydroxy group of a polyprenol with 1 mol eq. of phosphoric acid.

This implementation uses simplified heuristics:
  – It checks that the molecule contains a phosphate (or diphosphate) group.
  – It then looks for an “allylic phosphate” linkage where the phosphate is connected 
    via an oxygen to a carbon that, while not directly bearing a double bond, is adjacent 
    to a carbon that is involved in a C=C double bond.
  – It also requires that the molecule has a sufficiently long carbon chain with multiple 
    C=C bonds, as expected for a polyprenol.
If these conditions are met, the function returns True with an explanation.
Otherwise, it returns False with a reason.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polyprenol_phosphate(smiles: str):
    """
    Determines if a molecule belongs to the polyprenol phosphate class.
    A polyprenol phosphate is defined as a prenol phosphate arising from the condensation
    of the terminal allylic hydroxy group of a polyprenol with phosphoric acid.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is a polyprenol phosphate, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Heuristic 1: Look for an allylic phosphate (or diphosphate) linkage.
    # For each phosphorus (atomic number 15) in the molecule,
    # inspect its oxygen neighbors. For any such oxygen, check whether it is attached to a carbon.
    # That carbon (the candidate) must have at least one other carbon neighbor (besides the oxygen)
    # that is involved in a double bond (with some other carbon) which makes the candidate an allylic center.
    found_allylic_linkage = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 15:  # phosphorus atom (P)
            for o_atom in atom.GetNeighbors():
                if o_atom.GetAtomicNum() != 8:  # must be oxygen
                    continue
                # Get neighbors of the oxygen (expecting at least one carbon)
                for cand in o_atom.GetNeighbors():
                    if cand.GetIdx() == atom.GetIdx():
                        continue  # skip going back to phosphorus
                    if cand.GetAtomicNum() != 6:  # candidate must be carbon
                        continue
                    # Now, check if candidate carbon has at least one OTHER carbon neighbor that is
                    # involved in a double bond.
                    for neighbor in cand.GetNeighbors():
                        if neighbor.GetIdx() == o_atom.GetIdx():
                            continue  # skip the oxygen we came from
                        if neighbor.GetAtomicNum() != 6:
                            continue  # only interested in carbon neighbors
                        # Instead of checking the direct bond between cand and neighbor (which is single
                        # in an allylic situation), we inspect the neighbor’s bonds: if neighbor is involved
                        # in a double bond to some other carbon (other than cand), then cand is allylic.
                        for bond in neighbor.GetBonds():
                            # Get the other atom in the bond
                            begin_idx = bond.GetBeginAtomIdx()
                            end_idx = bond.GetEndAtomIdx()
                            other_idx = begin_idx if end_idx == neighbor.GetIdx() else end_idx
                            if other_idx == cand.GetIdx():
                                continue
                            other_atom = mol.GetAtomWithIdx(other_idx)
                            if other_atom.GetAtomicNum() != 6:
                                continue
                            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                                found_allylic_linkage = True
                                break
                        if found_allylic_linkage:
                            break
                    if found_allylic_linkage:
                        break
                if found_allylic_linkage:
                    break
        if found_allylic_linkage:
            break

    if not found_allylic_linkage:
        return False, "No allylic phosphate linkage found (phosphate not connected via an oxygen to an allylic carbon)"

    # Heuristic 2: Check that the molecule has a sufficiently long polyprenol chain.
    # Count the number of carbon atoms.
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 10:
        return False, f"Too few carbon atoms ({c_count}) for a polyprenol chain."
    
    # Heuristic 3: Count the number of C=C double bonds (signature of isoprene units).
    double_bond_count = 0
    for bond in mol.GetBonds():
        if (bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and 
            bond.GetBeginAtom().GetAtomicNum() == 6 and 
            bond.GetEndAtom().GetAtomicNum() == 6):
            double_bond_count += 1
    if double_bond_count < 2:
        return False, f"Not enough C=C bonds ({double_bond_count}) to support a polyprenol structure."
    
    # Optional: Check that the molecular weight is within typical range.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, f"Molecular weight ({mol_wt:.1f} Da) seems too low for a polyprenol phosphate."

    return True, "Contains polyprenol chain with multiple C=C bonds and an allylic phosphate (or diphosphate) linkage."

# Example usage (for testing):
if __name__ == '__main__':
    # Using one of the provided test SMILES as an example:
    test_smiles = "CC(C)=CCC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C/CC\\C(C)=C/CC\\C(C)=C/CC\\C(C)=C/CC\\C(C)=C/CC\\C(C)=C/CC\\C(C)=C/COP(O)(O)=O"
    result, reason = is_polyprenol_phosphate(test_smiles)
    print("Result:", result)
    print("Reason:", reason)