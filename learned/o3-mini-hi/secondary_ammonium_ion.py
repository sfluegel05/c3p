"""
Classifies: CHEBI:137419 secondary ammonium ion
"""
"""
Classifies: Secondary ammonium ion
Definition:
  "An organic cation obtained by protonation of any secondary amino compound; major species at pH 7.3."

A protonated secondary amine (secondary ammonium ion) arises from a secondary amine (R2NH)
by adding an extra proton to give R2NH2+. In this state the nitrogen:
  - Carries a formal charge of +1.
  - Has 4 bonds in total.
  - Has exactly 2 bonds to heavy (non-hydrogen) atoms and 2 bonds to hydrogen atoms.
  - Is sp3-hybridized.
  - Is not directly attached (through one of its heavy neighbors that is a carbon) 
    to a carbonyl group (a double-bonded oxygen).
Additionally, the overall molecule must be a cation (net charge > 0).
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_secondary_ammonium_ion(smiles: str):
    """
    Determines if the molecule (given by a SMILES string) contains at least one candidate secondary ammonium ion center.
    
    A candidate nitrogen atom must satisfy:
    - Be a nitrogen (atomic number 7) with a formal charge of +1.
    - Have a total of 4 bonds (after adding explicit hydrogens).
    - Have exactly 2 heavy neighbors (neighbors with atomic number != 1).
    - Have exactly 2 hydrogen atoms (counted via GetTotalNumHs()).
    - Be sp3 hybridized.
    - Not be attached (via any heavy neighbor that is carbon) 
      to a carbon that is double bonded to oxygen (i.e. a carbonyl group).
    - The overall molecule must be a cation (net charge > 0).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if a candidate secondary ammonium ion center is found, else False.
        str: Explanation for the decision.
    """
    # Parse the molecule from the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so we can count hydrogens properly.
    mol = Chem.AddHs(mol)
    
    # Check that the overall molecule is a cation (net formal charge > 0)
    net_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    if net_charge <= 0:
        return False, "Overall molecule is not a cation (net charge <= 0)."
    
    # Loop over every atom to find any candidate secondary ammonium ion center.
    for atom in mol.GetAtoms():
        # Check that the atom is nitrogen with exactly +1 formal charge.
        if atom.GetAtomicNum() != 7 or atom.GetFormalCharge() != 1:
            continue
        
        # Check that the nitrogen is sp3 hybridized.
        if atom.GetHybridization() != rdchem.HybridizationType.SP3:
            continue

        # After adding explicit hydrogens, degree reflects all bonds.
        total_neighbors = atom.GetDegree()
        if total_neighbors != 4:
            continue
        
        # Count the total number of hydrogen atoms (explicitly via GetTotalNumHs())
        hydrogen_count = atom.GetTotalNumHs()
        if hydrogen_count != 2:
            continue
        
        # Count heavy (non-hydrogen) neighbors.
        heavy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() != 1]
        if len(heavy_neighbors) != 2:
            continue
        
        # For each heavy neighbor that is carbon, check that it is not bound via a double bond to oxygen.
        violates_carbonyl_rule = False
        for nbr in heavy_neighbors:
            if nbr.GetAtomicNum() != 6:
                continue
            # Loop over bonds in the neighbor.
            for bond in nbr.GetBonds():
                # Skip the bond to our candidate nitrogen.
                if bond.GetOtherAtom(nbr).GetIdx() == atom.GetIdx():
                    continue
                # Check if the bond is a double bond
                if bond.GetBondType() == rdchem.BondType.DOUBLE:
                    other_atom = bond.GetOtherAtom(nbr)
                    if other_atom.GetAtomicNum() == 8:
                        # Found a carbonyl attachment.
                        violates_carbonyl_rule = True
                        break
            if violates_carbonyl_rule:
                break
        
        if violates_carbonyl_rule:
            continue  # Disqualify this candidate
        
        # If we get here, we found a candidate nitrogen that meets all criteria.
        return True, "Found a candidate secondary ammonium ion center matching criteria."
    
    # If no candidate nitrogen qualifies.
    return False, "No secondary ammonium ion center found based on the defined criteria."

# When run as a script, perform a few test cases.
if __name__ == "__main__":
    test_smiles_list = [
        "[H]\\C(C)=C1/C[NH2+][C@@]2([H])CC3=C(NC4=CC=CC=C34)C(=O)C[C@]1([H])[C@]2([H])C(=O)OC",  # perivine(1+)
        "C1COCC[NH2+]1",  # morpholinium
        "CC(C[NH2+]C)C",  # N-methyl-2-methylpropan-1-aminium
        "[C@]12([C@]3([C@@]([C@@]4(C(C[C@@H](O)CC4)=CC3)C)(CC[C@]2(C)[C@]5([C@@H]([C@]6(O[C@]5(C1)[H])CC[C@@H](C)C[NH2+]6)C)[H])[H])[H])[H]",  # solasodine(1+)
        "C1CCCC[NH2+]1",  # piperidinium
        "C[NH2+]C",  # dimethylaminium
        "C(CC[NH3+])CC[NH+](C)C",  # N,N-dimethyl-cadaverine(+2) -- candidate secondary site should be found here
    ]
    
    for s in test_smiles_list:
        result, reason = is_secondary_ammonium_ion(s)
        print("SMILES: {}\nResult: {}\nReason: {}\n".format(s, result, reason))