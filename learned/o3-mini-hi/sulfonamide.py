"""
Classifies: CHEBI:35358 sulfonamide
"""
#!/usr/bin/env python
"""
Classifies: Sulfonamide 
Definition: An amide of a sulfonic acid, RS(=O)2NR'2.
Criteria enforced here:
  – The sulfur atom (S) must be tetravalent (with exactly 4 bonds).
  – It must have exactly two oxygen neighbors connected via double bonds.
  – It must have exactly two other neighbors: one must be a carbon and the other an acyclic nitrogen.
This procedure tries to avoid false matching where the S(=O)(=O)–N motif appears in extended or inappropriate connectivity.
"""

from rdkit import Chem

def is_sulfonamide(smiles: str):
    """
    Determines if a molecule is a sulfonamide based on its SMILES string.
    A sulfonamide is defined as RS(=O)(=O)NR'2. In our approach, we look for a sulfur atom that:
      - Has exactly 4 connections.
      - Has exactly 2 bonds to oxygen with bond order 2.
      - Is also connected by single bonds to one carbon (R group) and one nitrogen.
      - The nitrogen must not be in a ring.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule contains a valid sulfonamide group, False otherwise.
        str: Explanation for the classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Iterate over atoms looking for a valid S (sulfur) candidate.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 16:  # We want sulfur only.
            continue
        # For a typical sulfonamide sulfur, the degree (number of bonds) should be 4.
        if atom.GetDegree() != 4:
            continue
        
        # Counters and flags:
        double_oxy_count = 0  # Count oxygen atoms double-bonded to S.
        carbon_found = False
        nitrogen_found = False
        n_neighbor = None
        
        # Check each bond/neighbour.
        for bond in atom.GetBonds():
            nb = bond.GetOtherAtom(atom)
            # Check bond order (we compare bond.GetBondTypeAsDouble() against 2.0 for a double bond)
            bond_order = bond.GetBondTypeAsDouble()
            if nb.GetAtomicNum() == 8 and abs(bond_order - 2.0) < 1e-3:
                double_oxy_count += 1
            # For the remaining substituents we require one nitrogen and one carbon:
            elif nb.GetAtomicNum() == 7 and bond_order == 1:
                # Check that this nitrogen is not part of a ring.
                if not nb.IsInRing():
                    nitrogen_found = True
                    n_neighbor = nb
            elif nb.GetAtomicNum() == 6 and bond_order == 1:
                carbon_found = True
            # If there are additional neighbors (e.g. another heteroatom) we do not want to count this S.
        
        # Now check that the sulfur atom meets our overall criteria.
        if double_oxy_count == 2 and nitrogen_found and carbon_found:
            # We also want to be sure that the S does not have any extra substituents.
            # In a perfect sulfonamide, S should have exactly 4 bonds: 2 double-bonded O, 1 C and 1 N.
            # (We already ensured degree==4.)
            return True, ("Molecule contains a sulfonamide group: RS(=O)(=O)-N (sulfur with two double-bonded oxygens, "
                          "bonded to a carbon and an acyclic nitrogen)")
    
    # If no such sulfur atom is found then the molecule is not a sulfonamide.
    return False, ("Molecule does not contain any sulfonamide substructure with the proper RS(=O)(=O)-N connectivity.")

# Simple testing if running as a script.
if __name__ == "__main__":
    # Test with one of the true positive examples from the outcomes (dofetilide)
    test_smiles = "CN(CCOc1ccc(NS(C)(=O)=O)cc1)CCc1ccc(NS(C)(=O)=O)cc1"  
    result, reason = is_sulfonamide(test_smiles)
    print("Result:", result)
    print("Reason:", reason)