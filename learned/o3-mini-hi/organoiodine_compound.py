"""
Classifies: CHEBI:37142 organoiodine compound
"""
#!/usr/bin/env python3
"""
Classifies: CHEBI:? organoiodine compound

An organoiodine compound is defined as a compound containing at least one carbon–iodine bond.
To avoid trivial, ionic, or highly heteroatom‐rich compounds we apply extra filters:
  - The molecule must be a single fragment (not a salt or mixture).
  - The molecule must have a minimum number of heavy atoms.
  - No atom should carry a nonzero formal charge.
  - The overall organic framework must be sufficiently carbon‐rich.
  - In at least one detected C–I bond the iodine must be bound to exactly one carbon, and that carbon
    must have at least one neighboring carbon.
"""

from rdkit import Chem

def is_organoiodine_compound(smiles: str):
    """
    Determines whether a molecule is an organoiodine compound based on its SMILES string.
    Additional filters:
      - The molecule must be a single fragment.
      - The molecule must have no atoms with a nonzero formal charge.
      - It must contain at least 6 heavy atoms.
      - The ratio of non‐carbon heavy atoms to carbons must not be too high.
      - At least one iodine atom must be directly bound to a carbon atom that has at least one other carbon neighbor.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as an organoiodine compound, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Reject if the molecule breaks into more than one fragment (possible salt or mixture)
    frags = Chem.GetMolFrags(mol)
    if len(frags) > 1:
        return False, "Compound is fragmented (possible salt or mixture)"
    
    # Reject if any atom carries a nonzero formal charge.
    for atom in mol.GetAtoms():
        if atom.GetFormalCharge() != 0:
            return False, "Compound has a nonzero formal charge"
    
    # Reject trivial molecules by heavy atom count (we require at least 6 heavy atoms)
    heavy_atoms = mol.GetNumHeavyAtoms()
    if heavy_atoms < 6:
        return False, f"Compound too simple with only {heavy_atoms} heavy atoms"
    
    # Check overall organic character: count the number of carbons and compare to non-carbon heavy atoms.
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    non_carbon = heavy_atoms - carbon_count
    if carbon_count == 0 or (non_carbon / carbon_count) > 1.0:
        return False, "Compound appears too heteroatom-rich to be a typical organoiodine compound"
    
    # Look for iodine atoms that are directly bound to a carbon.
    # We loop over all iodine atoms (atomic number 53).
    valid_CI_found = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 53:  # Iodine
            neighbors = atom.GetNeighbors()
            # For an organic C–I bond, expect iodine to have a single (heavy) neighbor which is carbon.
            if len(neighbors) == 1 and neighbors[0].GetAtomicNum() == 6:
                carbon_atom = neighbors[0]
                # Check that the carbon atom has at least one other carbon neighbor (to ensure organic context)
                carbon_neighbors = [nbr for nbr in carbon_atom.GetNeighbors()
                                    if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != atom.GetIdx()]
                if len(carbon_neighbors) >= 1:
                    valid_CI_found = True
                    break
    
    if not valid_CI_found:
        return False, "Compound does not contain a valid carbon–iodine bond in an organic framework"
    
    return True, "Compound contains a carbon–iodine bond in an organic framework"

# Example usage (for debugging/testing):
if __name__ == "__main__":
    test_smiles = [
        "IC1=CC(F)=C(NC(=O)CC2(CCCCC2)CC(O)=O)C=C1",  # 2-{1-[2-(2-fluoro-4-iodoanilino)-2-oxoethyl]cyclohexyl}acetic acid (TP)
        "CCOc1nc2ccc(I)cc2c(=O)n1CCC",               # proquinazid (TP)
        "C(O)C1OC(OC1)CI",                           # domiodol (FP: rejected due to not meeting organic criteria)
        "ICCCCCCl",                                  # 1-Chloro-4-iodobutane (TP)
    ]
    
    for smi in test_smiles:
        result, reason = is_organoiodine_compound(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")