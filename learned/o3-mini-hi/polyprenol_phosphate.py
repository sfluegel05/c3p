"""
Classifies: CHEBI:16460 polyprenol phosphate
"""
"""
Classifies: Polyprenol phosphate
Definition: A prenol phosphate resulting from the formal condensation of the terminal allylic hydroxy group of a polyprenol with 1 mol eq. of phosphoric acid.
This program uses simplified heuristics:
  – It parses the molecule.
  – It checks for the presence of a phosphorus atom (or diphosphate group)
    that is connected via an oxygen to an allylic carbon (a carbon that has at least one double bond to another carbon).
  – It also requires that the molecule has a significant carbon chain and contains at least a couple of C=C double bonds.
If all these conditions are met, it returns True along with an “explanation”.
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
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Heuristic 1: Check for an allylic phosphate linkage.
    # Look for phosphorus atoms (atomic number 15) in the molecule.
    found_allylic_linkage = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 15:  # phosphorus found
            # iterate over its neighbors to find oxygen atoms (which form the phosphate)
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 8:  # oxygen atom
                    # check if this oxygen is also bonded to a carbon (the linkage to the polyprenol)
                    for second_nbr in nbr.GetNeighbors():
                        # avoid going back to phosphorus
                        if second_nbr.GetIdx() == atom.GetIdx():
                            continue
                        # If the neighbor is a carbon, check if it is allylic.
                        if second_nbr.GetAtomicNum() == 6:
                            # Check if the carbon has at least one double bond to another carbon.
                            for bond in second_nbr.GetBonds():
                                # Get the other atom in the bond
                                other = bond.GetOtherAtom(second_nbr)
                                if other.GetAtomicNum() == 6 and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
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
    
    # Heuristic 2: Check that the molecule has a sufficiently long carbon chain.
    # We do a simple global count of carbon atoms.
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 10:
        return False, f"Too few carbon atoms ({c_count}) for a polyprenol chain."
    
    # Heuristic 3: Count the number of C=C double bonds (typical in isoprene units).
    double_bond_count = 0
    for bond in mol.GetBonds():
        if (bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and 
            bond.GetBeginAtom().GetAtomicNum() == 6 and 
            bond.GetEndAtom().GetAtomicNum() == 6):
            double_bond_count += 1
    if double_bond_count < 2:
        return False, f"Not enough C=C bonds ({double_bond_count}) to support a polyprenol structure."
        
    # Optionally, one might check the molecular weight to filter out very small molecules.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    # Note: polyprenols (and their phosphate derivatives) are often >300 Da; adjust as needed.
    if mol_wt < 300:
        return False, f"Molecular weight ({mol_wt:.1f} Da) seems too low for a polyprenol phosphate."
    
    # If all conditions are met, classify the structure as a polyprenol phosphate.
    return True, "Contains a long polyprenol chain with multiple C=C bonds and an allylic phosphate linkage."
    
# Example usage (for testing purpose):
if __name__ == '__main__':
    # One example given in the prompt: tritrans,heptacis-undecaprenyl phosphate
    test_smiles = "CC(C)=CCC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C/CC\\C(C)=C/CC\\C(C)=C/CC\\C(C)=C/CC\\C(C)=C/CC\\C(C)=C/CC\\C(C)=C/COP(O)(O)=O"
    result, reason = is_polyprenol_phosphate(test_smiles)
    print("Result:", result)
    print("Reason:", reason)