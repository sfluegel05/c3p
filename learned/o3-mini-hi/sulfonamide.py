"""
Classifies: CHEBI:35358 sulfonamide
"""
"""
Classifies: Sulfonamide
Definition: A sulfonamide is an amide of a sulfonic acid, i.e. its functional group has the form
            RS(=O)(=O)-N(R) (with the sulfur having exactly 4 neighbors: 2 double bonds to oxygen,
            one single (or aromatic) bond to nitrogen, and one single (or aromatic) bond to another atom).
We improved on the previous method by:
  • Using explicit bond type comparisons (checking both SINGLE and AROMATIC as “single‐bond equivalents”)
  • Not relying on floating point tolerances for bond orders.
  • Rejecting any S center that does not have exactly 4 neighbors.
This approach still may miss some edge cases but improves the balance between false positives and false negatives.
"""

from rdkit import Chem

def is_sulfonamide(smiles: str):
    """
    Determines if a molecule is a sulfonamide based on its SMILES string.
    The method looks for a sulfur (S) atom that is bonded to exactly:
      • Two oxygen atoms via a double bond (BondType.DOUBLE)
      • One nitrogen atom via a single or aromatic bond (BondType.SINGLE or BondType.AROMATIC)
      • One additional non-nitrogen heavy atom via a single or aromatic bond.
      
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if at least one sulfonamide group (RS(=O)(=O)-N(R)) is found, False otherwise.
        str: Explanation for the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Iterate over atoms; look for sulfur atoms.
    for atom in mol.GetAtoms():
        # Consider only sulfur atoms (atomic number 16)
        if atom.GetAtomicNum() != 16:
            continue
        
        # The sulfur in a sulfonamide should be tetra-coordinated.
        # (Note: GetDegree returns the number of explicit (bonded) neighbors.)
        if atom.GetDegree() != 4:
            continue
        
        # Initialize counters for the substituents:
        oxygen_double_bond = 0     # Count bonds that are double and to oxygen.
        nitrogen_single_bond = 0     # Count bonds that are single (or aromatic) to nitrogen.
        other_single_bond = 0        # Count bonds that are single (or aromatic) to non-nitrogen.
        valid = True  # flag to mark if any unexpected bond type is encountered.
        
        # Loop over all bonds for this S atom.
        for bond in atom.GetBonds():
            nbr = bond.GetOtherAtom(atom)
            bond_type = bond.GetBondType()
            
            # If bond is double and neighbor is oxygen then count.
            if bond_type == Chem.rdchem.BondType.DOUBLE and nbr.GetAtomicNum() == 8:
                oxygen_double_bond += 1
            # If bond is single or aromatic then count accordingly.
            elif bond_type in (Chem.rdchem.BondType.SINGLE, Chem.rdchem.BondType.AROMATIC):
                if nbr.GetAtomicNum() == 7:
                    nitrogen_single_bond += 1
                else:
                    other_single_bond += 1
            else:
                # If any bond is not a double, single, or aromatic bond, skip this S atom.
                valid = False
                break
        
        # Verify that this sulfur center has exactly the following connectivity:
        # 2 double bonds to oxygen, 1 bond to nitrogen, 1 bond to non-nitrogen.
        if valid and oxygen_double_bond == 2 and nitrogen_single_bond == 1 and other_single_bond == 1:
            return True, ("Molecule contains a sulfonamide group: RS(=O)(=O)-N(R) "
                          "(sulfur with exactly two double-bonded oxygens, one single/aromatic bond to nitrogen, "
                          "and one single/aromatic bond to another substituent)")
    
    # If no matching S center was found, return False.
    return False, ("Molecule does not contain a sulfonamide substructure with the required "
                   "RS(=O)(=O)-N connectivity.")

# For simple testing when running as a script.
if __name__ == "__main__":
    # Test with an example known sulfonamide (e.g., dofetilide).
    test_smiles = "CN(CCOc1ccc(NS(C)(=O)=O)cc1)CCc1ccc(NS(C)(=O)=O)cc1"
    result, reason = is_sulfonamide(test_smiles)
    print("Result:", result)
    print("Reason:", reason)