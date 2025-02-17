"""
Classifies: CHEBI:35358 sulfonamide
"""
"""
Classifies: Sulfonamide
Definition: A sulfonamide is an amide of a sulfonic acid, i.e. it has the form
            RS(=O)₂NR'₂.
In our implementation we require that a candidate sulfonamide group satisfy:
  • The sulfur (atomic number 16) is tetravalent (degree==4).
  • Exactly two of its bonds are double bonds to oxygen atoms.
  • Its remaining two bonds (both single bonds) go to one carbon and one nitrogen.
Note: Although many texts require the sulfonamide nitrogen to be acyclic, our
      analysis of the outcomes shows that excluding cyclic nitrogens causes
      false negatives. (Some valid sulfonamides contain a cyclic substructure.)
      
An accepted match will return True with a corresponding explanation; if no
such S atom is found the function returns False.
"""

from rdkit import Chem

def is_sulfonamide(smiles: str):
    """
    Determines if a molecule is a sulfonamide based on its SMILES string.
    A sulfonamide is defined as RS(=O)(=O)NR'₂. In this approach we iterate over
    all sulfur (S) atoms in the molecule and check the following:
      - The sulfur atom has exactly 4 bonds.
      - Exactly two of these bonds are double bonds to oxygen (O) atoms.
      - Exactly one remaining bond is a single bond to a carbon (C) atom.
      - Exactly one remaining bond is a single bond to a nitrogen (N) atom.
    (Any additional substituents result in rejecting that particular S atom.)
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if a sulfonamide group is present, False otherwise.
        str: Explanation for the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Iterate over atoms and search for a valid sulfonamide S candidate.
    for atom in mol.GetAtoms():
        # Look only at sulfur atoms.
        if atom.GetAtomicNum() != 16:
            continue
        # Check that the sulfur is tetravalent (degree==4).
        if atom.GetDegree() != 4:
            continue
        
        # Set counters for the substituents.
        double_bonded_oxygen = 0
        single_bonded_carbon = 0
        single_bonded_nitrogen = 0
        
        # Iterate over bonds from this sulfur atom.
        for bond in atom.GetBonds():
            nbr = bond.GetOtherAtom(atom)
            bond_order = bond.GetBondTypeAsDouble()
            if nbr.GetAtomicNum() == 8 and abs(bond_order - 2.0) < 1e-3:
                double_bonded_oxygen += 1
            elif nbr.GetAtomicNum() == 6 and abs(bond_order - 1.0) < 1e-3:
                single_bonded_carbon += 1
            elif nbr.GetAtomicNum() == 7 and abs(bond_order - 1.0) < 1e-3:
                single_bonded_nitrogen += 1
            else:
                # If the neighbor is not O (double bond), C (single bond) or N (single bond),
                # then this sulfur is not the expected sulfonamide center.
                pass
        
        # Now check if the connectivity exactly matches:
        # Two double-bonded oxygens, one single-bonded carbon and one single-bonded nitrogen.
        if double_bonded_oxygen == 2 and single_bonded_carbon == 1 and single_bonded_nitrogen == 1:
            return True, ("Molecule contains a sulfonamide group: RS(=O)(=O)-N(R) "
                          "(sulfur with two double-bonded oxygens bonded to one carbon and one nitrogen)")
    
    # If no matching sulfonamide group is found return False.
    return False, ("Molecule does not contain a sulfonamide substructure with the proper RS(=O)(=O)-N connectivity.")

# The following can be used for simple testing when running as a script.
if __name__ == "__main__":
    # For example, test with dofetilide SMILES.
    test_smiles = "CN(CCOc1ccc(NS(C)(=O)=O)cc1)CCc1ccc(NS(C)(=O)=O)cc1"  
    result, reason = is_sulfonamide(test_smiles)
    print("Result:", result)
    print("Reason:", reason)