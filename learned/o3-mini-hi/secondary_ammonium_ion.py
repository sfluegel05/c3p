"""
Classifies: CHEBI:137419 secondary ammonium ion
"""
"""
Classifies: Secondary ammonium ion
Definition:
   "An organic cation obtained by protonation of any secondary amino compound;
    major species at pH 7.3."

A protonated secondary amine (secondary ammonium ion) arises from a secondary amine (R2NH)
by adding an extra proton to give R2NH2+. In this state the nitrogen:
  - Carries a formal charge of +1.
  - Has 4 bonds in total.
  - Has exactly 2 bonds to heavy atoms (non-hydrogen atoms) and 2 bonds to hydrogens.
  - Is not directly attached (through a heavy neighbor) to an “amidic” group 
    (i.e. a heavy neighbor that has a double bond to oxygen on another bond).
Additionally, we require that the overall molecule is a cation (net charge > 0).
"""

from rdkit import Chem

def is_secondary_ammonium_ion(smiles: str):
    """
    Determines if the molecule (given by a SMILES string) contains at least one secondary ammonium ion center.
    
    For a candidate nitrogen atom to qualify as a secondary ammonium ion center, it must:
      - Be nitrogen (atomic number 7) with a formal charge of +1.
      - After adding explicit hydrogens, have exactly two heavy (non-hydrogen) neighbors.
      - Have exactly two hydrogens attached.
      - Not be linked (via its heavy neighbors) to a carbonyl group—specifically, any heavy neighbor
        that is a carbon double bonded to oxygen (on a bond other than the one to the candidate nitrogen)
        disqualifies that candidate.
      - The overall molecule must have a net charge > 0.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if a candidate secondary ammonium ion center is found, False otherwise.
        str: Reason message detailing the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that we can accurately assess hydrogen counts.
    mol = Chem.AddHs(mol)
    
    # Ensure that the overall molecule is a cation (net formal charge > 0).
    net_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    if net_charge <= 0:
        return False, "Molecule overall is not a cation (net charge <= 0)."
    
    # Loop over atoms to check for a candidate secondary ammonium ion center.
    for atom in mol.GetAtoms():
        # Candidate must be nitrogen with a +1 formal charge.
        if atom.GetAtomicNum() != 7 or atom.GetFormalCharge() != 1:
            continue
        
        # Instead of relying on hybridization and degree, check the connectivity:
        # Count heavy (non-H) neighbors and hydrogens.
        heavy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() != 1]
        if len(heavy_neighbors) != 2:
            continue  # Should be bonded to exactly 2 heavy atoms.
        
        # Count total hydrogens attached (explicit hydrogens).
        total_h = atom.GetTotalNumHs()
        if total_h != 2:
            continue  # Should be bonded to exactly 2 hydrogens.
        
        # Screen out candidates if any heavy neighbor appears "amidic".
        # For each heavy neighbor, look at bonds (excluding the bond to the candidate nitrogen).
        amidic_flag = False
        for nbr in heavy_neighbors:
            # Iterate over bonds of the neighbor.
            for bond in nbr.GetBonds():
                # Skip the bond that connects this neighbor to our candidate nitrogen.
                if bond.GetOtherAtom(nbr).GetIdx() == atom.GetIdx():
                    continue
                # If the heavy neighbor is carbon and is double bonded to oxygen, mark as amidic.
                if nbr.GetAtomicNum() == 6:
                    other = bond.GetOtherAtom(nbr)
                    if other.GetAtomicNum() == 8 and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        amidic_flag = True
                        break
            if amidic_flag:
                break
        if amidic_flag:
            continue
        
        # If a candidate passes all these criteria, return True.
        return True, "Found a candidate secondary ammonium ion center matching the criteria."
    
    # If no candidate was found, return False with the appropriate reason.
    return False, "No secondary ammonium ion center found based on the defined criteria."

# Example test runs (run only if called as a standalone script).
if __name__ == "__main__":
    test_smiles = [
        "[H]\\C(C)=C1/C[NH2+][C@@]2([H])CC3=C(NC4=CC=CC=C34)C(=O)C[C@]1([H])[C@]2([H])C(=O)OC",  # perivine(1+)
        "C1COCC[NH2+]1",  # morpholinium
        "CC(C[NH2+]C)C",  # N-methyl-2-methylpropan-1-aminium
        "[C@]12([C@]3([C@@]([C@@]4(C(C[C@@H](O)CC4)=CC3)C)(CC[C@]2(C)[C@]5([C@@H]([C@]6(O[C@]5(C1)[H])CC[C@@H](C)C[NH2+]6)C)[H])[H])[H])[H]",  # solasodine(1+)
        "C1CCCC[NH2+]1",  # piperidinium
    ]
    
    for s in test_smiles:
        result, reason = is_secondary_ammonium_ion(s)
        print("SMILES: {}\nResult: {}\nReason: {}\n".format(s, result, reason))