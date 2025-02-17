"""
Classifies: CHEBI:27933 beta-lactam antibiotic
"""
"""
Classifies: Beta-lactam antibiotics
Definition: An organonitrogen heterocyclic antibiotic that contains a beta-lactam ring.
A beta-lactam ring is a four-membered cyclic amide that contains one nitrogen and three carbons,
with one of the carbons forming a carbonyl group (C=O).
"""

from rdkit import Chem

def is_beta_lactam_antibiotic(smiles: str):
    """
    Determines if a molecule is a beta-lactam antibiotic based on its SMILES string.
    We classify a molecule as a beta-lactam antibiotic if it contains at least one beta-lactam ring,
    which is defined as a 4-membered ring with one nitrogen, three carbons, and where one of the carbons
    is double-bonded to an oxygen (i.e. a carbonyl group).

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule contains a beta-lactam ring, False otherwise.
        str: Explanation for the classification.
    """
    # Parse the SMILES string to an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Retrieve all rings in the molecule.
    ring_info = mol.GetRingInfo().AtomRings()
    
    # Iterate over each ring.
    for ring in ring_info:
        # We are interested in 4-membered rings.
        if len(ring) == 4:
            # Count nitrogen and carbon atoms in the ring.
            n_count = 0
            c_atoms = []  # Store indices of carbon atoms in the ring.
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() == 7:
                    n_count += 1
                elif atom.GetAtomicNum() == 6:
                    c_atoms.append(idx)
                else:
                    # If any other element appears in the ring, skip this ring.
                    n_count = -999
                    break
            # Check if the ring meets our criteria: exactly one N, three C.
            if n_count == 1 and len(c_atoms) == 3:
                # Now check if at least one of the carbon atoms in the ring is a carbonyl carbon,
                # i.e. has a double bond to an oxygen.
                carbonyl_found = False
                for c_idx in c_atoms:
                    carbon = mol.GetAtomWithIdx(c_idx)
                    # Iterate over the neighbors of the carbon.
                    for neighbor in carbon.GetNeighbors():
                        # check if neighbor is oxygen
                        if neighbor.GetAtomicNum() == 8:
                            bond = mol.GetBondBetweenAtoms(c_idx, neighbor.GetIdx())
                            # Check if the bond is a double bond.
                            if bond is not None and bond.GetBondType() == Chem.BondType.DOUBLE:
                                carbonyl_found = True
                                break
                    if carbonyl_found:
                        break
                if carbonyl_found:
                    return True, "Molecule contains a beta-lactam ring (4-membered ring with 1 N, 3 C, one C=O)."
    
    return False, "No beta-lactam ring found in the molecule"