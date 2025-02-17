"""
Classifies: CHEBI:16460 polyprenol phosphate
"""
"""
Classifies: Polyprenol phosphate
Definition: A prenol phosphate resulting from the formal condensation of 
the terminal allylic hydroxy group of a polyprenol with 1 mol eq. of phosphoric acid.
This implementation uses simplified heuristics:
  – It checks that the molecule contains a phosphate group (or diphosphate) linked via an oxygen.
  – Instead of requiring the carbon directly attached to that oxygen to be sp2 or have a double bond,
    it looks at the neighbors of that carbon to see if one of them is involved in a C=C double bond.
  – It also requires that the molecule has a sufficiently long carbon chain and multiple double bonds.
If these conditions are met, the function returns True with an explanation.
Otherwise, it returns False with a reason.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polyprenol_phosphate(smiles: str):
    """
    Determines if a molecule belongs to the polyprenol phosphate class.
    A polyprenol phosphate is defined as a prenol phosphate arising from the condensation
    of the terminal allylic hydroxy group of a polyprenol with 1 mol eq. of phosphoric acid.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is a polyprenol phosphate, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Heuristic 1: Look for an allylic phosphate linkage.
    # We search for a phosphorus atom (atomic number 15) attached to oxygen.
    # Then for each such O, we follow to the carbon it is bound to.
    # Finally, we check if that carbon is adjacent (via a single bond) to another carbon
    # that is involved in a carbon-carbon double bond.
    found_allylic_linkage = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 15:  # phosphorus atom
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 8:  # oxygen neighbor
                    # Now, from the oxygen, find a neighbor that is a carbon (and not the phosphorus)
                    for carbon in nbr.GetNeighbors():
                        if carbon.GetIdx() == atom.GetIdx():
                            continue
                        if carbon.GetAtomicNum() != 6:
                            continue
                        # Instead of requiring 'carbon' itself to have a double bond,
                        # check if any of its neighbors (except the oxygen) is involved in a C=C double bond.
                        for allyl_candidate in carbon.GetNeighbors():
                            if allyl_candidate.GetIdx() == nbr.GetIdx():
                                continue  # skip back to the oxygen
                            if allyl_candidate.GetAtomicNum() != 6:
                                continue
                            # Now, check the bond between carbon and allyl_candidate to see if it is a double bond.
                            bond = mol.GetBondBetweenAtoms(carbon.GetIdx(), allyl_candidate.GetIdx())
                            if bond and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                                found_allylic_linkage = True
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
    # We use a simple count of carbon atoms.
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 10:
        return False, f"Too few carbon atoms ({c_count}) for a polyprenol chain."
    
    # Heuristic 3: Count the number of C=C double bonds (a signature of isoprene units).
    double_bond_count = 0
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        if (bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and 
            a1.GetAtomicNum() == 6 and 
            a2.GetAtomicNum() == 6):
            double_bond_count += 1
    if double_bond_count < 2:
        return False, f"Not enough C=C bonds ({double_bond_count}) to support a polyprenol structure."
        
    # Optional molecular weight test.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, f"Molecular weight ({mol_wt:.1f} Da) seems too low for a polyprenol phosphate."
    
    # If all conditions are met, we classify the structure as a polyprenol phosphate.
    return True, "Contains a polyprenol chain with multiple C=C bonds and an allylic phosphate (or diphosphate) linkage."
    
# Example usage (for testing purposes):
if __name__ == '__main__':
    # Try one of the provided test SMILES: tritrans,heptacis-undecaprenyl phosphate
    test_smiles = "CC(C)=CCC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C/CC\\C(C)=C/CC\\C(C)=C/CC\\C(C)=C/CC\\C(C)=C/CC\\C(C)=C/CC\\C(C)=C/COP(O)(O)=O"
    result, reason = is_polyprenol_phosphate(test_smiles)
    print("Result:", result)
    print("Reason:", reason)