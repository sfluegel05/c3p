"""
Classifies: CHEBI:27933 beta-lactam antibiotic
"""
"""
Classifies: Beta-lactam antibiotics
Definition: An organonitrogen heterocyclic antibiotic that contains a beta-lactam ring.
A beta-lactam ring is defined as a four-membered cyclic amide containing exactly one nitrogen,
three carbons, and with one of the carbons bearing a double bond to oxygen.
This implementation looks for a 4-membered ring (using RDKit’s ring info) that satisfies these criteria.
"""

from rdkit import Chem

def is_beta_lactam_antibiotic(smiles: str):
    """
    Determines if a molecule appears to be a beta-lactam antibiotic based on its SMILES string.
    The molecule is considered to be a beta-lactam antibiotic if it contains a four-membered ring
    (cyclic amide) that:
      - consists of exactly one nitrogen and three carbons, and
      - has at least one carbon atom that is double-bonded to an oxygen (C=O)
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if a beta-lactam ring is found, False otherwise.
        str: Explanation for the classification.
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Attempt to standardize bond orders by kekulizing the molecule.
    try:
        Chem.Kekulize(mol, clearAromaticFlags=True)
    except Exception:
        # Proceed even if kekulization fails.
        pass

    # Get ring information (list of tuples, each tuple has atom indices forming a ring)
    ring_info = mol.GetRingInfo().AtomRings
    beta_lactam_found = False

    # Iterate over each ring
    for ring in ring_info:
        # Consider only rings with exactly 4 atoms.
        if len(ring) != 4:
            continue

        # Get the atom objects for the ring atoms.
        atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        # Count number of nitrogens and carbons in ring.
        n_count = sum(1 for atom in atoms if atom.GetSymbol() == 'N')
        c_count = sum(1 for atom in atoms if atom.GetSymbol() == 'C')
        
        # The beta-lactam ring must have exactly 1 N and 3 C.
        if n_count != 1 or c_count != 3:
            continue

        # Look for a carbonyl: at least one carbon in the ring must have a double bond to oxygen.
        carbonyl_found = False
        for atom in atoms:
            if atom.GetSymbol() != 'C':
                continue
            # Examine all bonds of this carbon atom.
            for bond in atom.GetBonds():
                # Check if the bond is a double bond.
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    # Check if the bonded neighbor is an oxygen.
                    neighbor = bond.GetOtherAtom(atom)
                    if neighbor.GetSymbol() == 'O':
                        carbonyl_found = True
                        break
            if carbonyl_found:
                break
        
        # If the ring passes both the atom count and carbonyl criteria, we assume a beta-lactam ring.
        if carbonyl_found:
            beta_lactam_found = True
            break

    if beta_lactam_found:
        return True, ("Molecule contains a beta-lactam ring (a 4-membered cyclic amide with 1 N, 3 C, "
                      "and one carbon bearing a C=O), which is typical for beta-lactam antibiotics.")
    else:
        return False, "No suitable beta-lactam ring found in the molecule"

# (Optional) Testing of some example SMILES strings.
if __name__ == "__main__":
    # Examples; note that a wide range of beta-lactam structures exist.
    example_smiles = [
        "[H][C@]12SCC(C)=C(N1C(=O)[C@H]2N)C(O)=O",  # 7beta-aminodeacetoxycephalosporanic acid (should be True)
        "O=C1CCN1",                                # azetidin-2-one (should be True)
        "C1=CC=CC=C1"                              # benzene (should be False)
    ]
    for smi in example_smiles:
        result, explanation = is_beta_lactam_antibiotic(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {explanation}\n")