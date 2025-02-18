"""
Classifies: CHEBI:35358 sulfonamide
"""
"""
Classifies: Sulfonamide
Definition: A sulfonamide is an amide of a sulfonic acid, i.e. it has the form
            RS(=O)₂NR'₂.
In our implementation we require that a candidate sulfonamide group satisfy:
  • The sulfur (atomic number 16) is tetravalent.
  • Exactly two of its bonds (with bond order near 2) are to oxygen atoms.
  • Exactly one of its single bonds (bond order near 1) is to a nitrogen atom.
  • The remaining single bond (bond order near 1) is to any other atom (typically an R group).
Any extra or unexpected bonds cause that sulfur to be rejected.
This approach does not attempt to evaluate the connectivity at the nitrogen.
"""

from rdkit import Chem

def is_sulfonamide(smiles: str):
    """
    Determines if a molecule is a sulfonamide based on its SMILES string.
    A sulfonamide is defined as RS(=O)(=O)NR'₂. In this approach we iterate over
    all sulfur (S) atoms in the molecule and check that at least one of them satisfies:
      - The S atom has exactly 4 bonds.
      - Exactly two bonds are double bonds (bond order near 2) to oxygen atoms.
      - Exactly one single bond (order near 1) goes to a nitrogen atom.
      - Exactly one single bond (order near 1) goes to a non-nitrogen atom.
      
    If such a sulfur is found, the function return True along with an explanation.
    Otherwise it returns False with an appropriate reason.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if a sulfonamide group is found, False otherwise.
        str: Explanation for the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Iterate over atoms; look for sulfur atoms.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 16:  # Only sulfur
            continue
        # Check that sulfur is tetravalent.
        if atom.GetDegree() != 4:
            continue
        
        # Counters for expected substituents:
        double_bonded_oxygen = 0
        single_bonded_nitrogen = 0
        other_single_bond = 0
        valid_candidate = True  # flag to mark any unexpected connectivity
        
        # Iterate over all bonds from the sulfur.
        for bond in atom.GetBonds():
            # Get the neighboring atom and bond order.
            nbr = bond.GetOtherAtom(atom)
            # Using GetBondTypeAsDouble helps to compare bond orders.
            bond_order = bond.GetBondTypeAsDouble()
            
            # Check if the bond is a double bond (order near 2) to an oxygen.
            if nbr.GetAtomicNum() == 8 and abs(bond_order - 2.0) < 1e-3:
                double_bonded_oxygen += 1
            # Check if the bond is a single bond (order near 1) to a nitrogen.
            elif nbr.GetAtomicNum() == 7 and abs(bond_order - 1.0) < 1e-3:
                single_bonded_nitrogen += 1
            # Check if the bond is a single bond (order near 1) to some other atom.
            elif abs(bond_order - 1.0) < 1e-3:
                other_single_bond += 1
            else:
                # Any other bond (e.g. bond orders that are not nearly 1 or 2)
                # or bonds to unexpected atoms cause this S to be rejected.
                valid_candidate = False
                break
        
        # Validate the candidate S atom:
        # We require exactly two double bonds to oxygen,
        # exactly one single bond to nitrogen,
        # and exactly one single bond to some other atom.
        if valid_candidate and double_bonded_oxygen == 2 and single_bonded_nitrogen == 1 and other_single_bond == 1:
            return True, ("Molecule contains a sulfonamide group: RS(=O)(=O)-N(R) "
                          "(sulfur with two double-bonded oxygens, one single-bonded nitrogen, "
                          "and one additional substituent)")
    
    # If no matching S center was found, return False with explanation.
    return False, ("Molecule does not contain a sulfonamide substructure with the required "
                   "RS(=O)(=O)-N connectivity.")

# For simple testing when running as a script.
if __name__ == "__main__":
    # Example: test with dofetilide SMILES which is known to be a sulfonamide.
    test_smiles = "CN(CCOc1ccc(NS(C)(=O)=O)cc1)CCc1ccc(NS(C)(=O)=O)cc1"
    result, reason = is_sulfonamide(test_smiles)
    print("Result:", result)
    print("Reason:", reason)