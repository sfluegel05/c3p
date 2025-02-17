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
  - Is not directly attached (through a heavy neighbor that is a carbon) to a carbonyl group (C(=O)).
Additionally, the overall molecule must be a cation (net charge > 0).
"""
from rdkit import Chem

def is_secondary_ammonium_ion(smiles: str):
    """
    Determines if the molecule (given by a SMILES string) contains at least one secondary ammonium ion center.
    
    For a candidate nitrogen atom to qualify as a secondary ammonium ion center, it must:
      - Be nitrogen (atomic number 7) with a formal charge of +1.
      - After adding explicit hydrogens, have exactly two heavy (non-hydrogen) neighbors.
      - Have exactly two hydrogen atoms attached (counted explicitly via neighbors).
      - Not be attached (through any heavy neighbor) to a carbon that is double bonded to oxygen.
      - The overall molecule must have a net charge > 0.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if a candidate secondary ammonium ion center is found, False otherwise.
        str: Reason message detailing the classification.
    """
    # Parse the molecule from the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that the number of hydrogen atoms can be counted directly.
    mol = Chem.AddHs(mol)
    
    # Check that the overall molecule is a cation (net formal charge > 0).
    net_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    if net_charge <= 0:
        return False, "Overall molecule is not a cation (net charge <= 0)."
    
    # Loop over all atoms to locate candidate secondary ammonium ion centers.
    for atom in mol.GetAtoms():
        # Candidate must be a nitrogen with a +1 formal charge.
        if atom.GetAtomicNum() != 7 or atom.GetFormalCharge() != 1:
            continue
        
        # Get heavy neighbors (non-hydrogen atoms).
        heavy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() != 1]
        # For a protonated secondary amine, expect exactly 2 heavy neighbors.
        if len(heavy_neighbors) != 2:
            continue
        
        # Count the number of hydrogen atoms attached directly.
        hydrogen_count = sum(1 for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 1)
        if hydrogen_count != 2:
            continue
        
        # Check each heavy neighbor to ensure it is not "amidic"
        # i.e. disqualify if any heavy neighbor (if it is carbon) has a double bond to oxygen
        # Note: we only check bonds other than the one to our candidate nitrogen.
        amidic_flag = False
        for nbr in heavy_neighbors:
            # Only consider carbon atoms.
            if nbr.GetAtomicNum() != 6:
                continue
            for bond in nbr.GetBonds():
                # Skip the bond that connects the heavy neighbor to our candidate.
                if bond.GetOtherAtom(nbr).GetIdx() == atom.GetIdx():
                    continue
                # Check if the bond is a double bond.
                if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    other = bond.GetOtherAtom(nbr)
                    # If the other atom in the bond is oxygen, mark as amidic.
                    if other.GetAtomicNum() == 8:
                        amidic_flag = True
                        break
            if amidic_flag:
                break
        if amidic_flag:
            continue
        
        # If a candidate nitrogen meets all the criteria, return True.
        return True, "Found a candidate secondary ammonium ion center matching criteria."
    
    # If no candidate nitrogen qualifies, return False with an explanation.
    return False, "No secondary ammonium ion center found based on the defined criteria."

# Example test runs (only executed when the script is run as main).
if __name__ == "__main__":
    test_smiles = [
        "[H]\\C(C)=C1/C[NH2+][C@@]2([H])CC3=C(NC4=CC=CC=C34)C(=O)C[C@]1([H])[C@]2([H])C(=O)OC",  # perivine(1+)
        "C1COCC[NH2+]1",  # morpholinium
        "CC(C[NH2+]C)C",  # N-methyl-2-methylpropan-1-aminium
        "[C@]12([C@]3([C@@]([C@@]4(C(C[C@@H](O)CC4)=CC3)C)(CC[C@]2(C)[C@]5([C@@H]([C@]6(O[C@]5(C1)[H])CC[C@@H](C)C[NH2+]6)C)[H])[H])[H])[H]",  # solasodine(1+)
        "C1CCCC[NH2+]1",  # piperidinium
        "C[NH2+]C",  # dimethylaminium
    ]
    
    for s in test_smiles:
        result, reason = is_secondary_ammonium_ion(s)
        print("SMILES: {}\nResult: {}\nReason: {}\n".format(s, result, reason))