"""
Classifies: CHEBI:27933 beta-lactam antibiotic
"""
"""
Classifies: Beta-lactam antibiotics
Definition: An organonitrogen heterocyclic antibiotic that contains a beta-lactam ring.
A beta-lactam ring is defined as a four-membered cyclic amide that contains exactly one nitrogen,
three carbons, and at least one carbon with a C=O double bond.
"""

from rdkit import Chem

def is_beta_lactam_antibiotic(smiles: str):
    """
    Determines if a molecule appears to be a beta-lactam antibiotic based on its SMILES string.
    The molecule is considered a beta-lactam antibiotic if it contains a four-membered ring
    (cyclic amide) that:
      - consists of exactly one nitrogen and three carbons, and
      - has at least one carbon atom in the ring that is double-bonded to an oxygen (C=O)
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if a valid beta-lactam ring is found, False otherwise.
        str: Explanation for the classification.
    """
    # Parse the SMILES string into a molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Attempt to kekulize - this may adjust bond orders for clearer matching.
    try:
        Chem.Kekulize(mol, clearAromaticFlags=True)
    except Exception:
        # Continue even if kekulization fails.
        pass

    # Retrieve the ring information by calling the AtomRings() method.
    ring_info = mol.GetRingInfo().AtomRings()
    beta_lactam_found = False

    # Iterate over each ring in the molecule.
    for ring in ring_info:
        # Only consider rings with exactly 4 atoms.
        if len(ring) != 4:
            continue
        
        # Retrieve the atom objects for the atoms of the ring.
        atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        
        # Count the number of nitrogen and carbon atoms in the ring.
        n_count = sum(1 for atom in atoms if atom.GetSymbol() == 'N')
        c_count = sum(1 for atom in atoms if atom.GetSymbol() == 'C')
        
        # A beta-lactam ring should have exactly 1 nitrogen and 3 carbons.
        if n_count != 1 or c_count != 3:
            continue

        # Look for a carbonyl group: at least one of the ring carbon atoms must have 
        # a double bond to an oxygen atom.
        carbonyl_found = False
        for atom in atoms:
            if atom.GetSymbol() != 'C':
                continue
            # Check each bond for a double bond to oxygen.
            for bond in atom.GetBonds():
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    neighbor = bond.GetOtherAtom(atom)
                    if neighbor.GetSymbol() == 'O':
                        carbonyl_found = True
                        break
            if carbonyl_found:
                break
        
        # If both criteria are met, we have a beta-lactam ring.
        if carbonyl_found:
            beta_lactam_found = True
            break

    if beta_lactam_found:
        return True, ("Molecule contains a beta-lactam ring (a 4-membered cyclic amide with 1 N, 3 C, "
                      "and a carbon with a C=O double bond), which is typical for beta-lactam antibiotics.")
    else:
        return False, "No suitable beta-lactam ring found in the molecule"

# (Optional) Test examples of beta-lactam antibiotics.
if __name__ == "__main__":
    test_smiles = [
        "[H][C@]12SCC(C)=C(N1C(=O)[C@H]2N)C(O)=O", # 7beta-aminodeacetoxycephalosporanic acid (should be True)
        "O=C1CCN1",                               # azetidin-2-one (should be True)
        "C1=CC=CC=C1"                             # benzene (should be False)
    ]
    for smi in test_smiles:
        result, explanation = is_beta_lactam_antibiotic(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {explanation}\n")