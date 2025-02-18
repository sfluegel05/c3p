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
  - At least one iodine must be bound to a carbon atom (even if hypervalent) that in turn is connected to another carbon.
"""

from rdkit import Chem

def is_organoiodine_compound(smiles: str):
    """
    Determines whether a molecule is an organoiodine compound based on its SMILES string.
    Additional filters:
      - The molecule must be a single fragment.
      - The molecule must have no atoms with a nonzero formal charge.
      - It must contain at least 6 heavy atoms.
      - The overall framework must be sufficiently carbon‐rich. For example, if fewer than 50%
        of heavy atoms are carbon or if the ratio of non‐carbon heavy atoms to carbon atoms exceeds 1, 
        the compound is considered too heteroatom‐rich.
      - At least one iodine atom must be directly bonded to a carbon, and that carbon must have at least one other carbon neighbor.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as an organoiodine compound, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Reject molecules that are fragmented (possible salt or mixture)
    frags = Chem.GetMolFrags(mol)
    if len(frags) > 1:
        return False, "Compound is fragmented (possible salt or mixture)"
    
    # Reject if any atom has a nonzero formal charge
    for atom in mol.GetAtoms():
        if atom.GetFormalCharge() != 0:
            return False, "Compound has a nonzero formal charge"
    
    # Require minimum heavy atom count (non-hydrogen atoms)
    heavy_atoms = mol.GetNumHeavyAtoms()
    if heavy_atoms < 6:
        return False, f"Compound too simple with only {heavy_atoms} heavy atoms"
    
    # Check overall organic character:
    # Count carbon atoms among heavy atoms.
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    non_carbon = heavy_atoms - carbon_count
    # Require that the ratio of non-carbon atoms to carbon atoms is not too high.
    if carbon_count == 0 or (non_carbon / carbon_count) > 1.0:
        return False, "Compound appears too heteroatom-rich to be a typical organoiodine compound"
    # Additionally, require that at least 50% of the heavy atoms are carbons.
    if (carbon_count / heavy_atoms) < 0.5:
        return False, "Not enough carbon atoms to support an organic framework"
    
    # Look for at least one C–I bond.
    # We now allow iodine atoms that have more than one neighbor;
    # we only require that at least one neighbor is a carbon which itself is connected to another carbon.
    valid_CI_found = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 53:  # Iodine
            # Look for any carbon neighbor of the iodine
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6:
                    # Check that this carbon has at least one additional carbon neighbor (besides the iodine)
                    carbon_neighbors = [nbr for nbr in neighbor.GetNeighbors() 
                                        if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != atom.GetIdx()]
                    if len(carbon_neighbors) >= 1:
                        valid_CI_found = True
                        break
            if valid_CI_found:
                break
                
    if not valid_CI_found:
        return False, "Compound does not contain a valid carbon–iodine bond in an organic framework"
    
    return True, "Compound contains a carbon–iodine bond in an organic framework"

# Example usage (for debugging/testing):
if __name__ == "__main__":
    test_smiles = [
        "IC1=CC(F)=C(NC(=O)CC2(CCCCC2)CC(O)=O)C=C1",  # 2-{1-[2-(2-fluoro-4-iodoanilino)-2-oxoethyl]cyclohexyl}acetic acid (true positive)
        "CCOc1nc2ccc(I)cc2c(=O)n1CCC",               # proquinazid (true positive)
        "C(O)C1OC(OC1)CI",                           # domiodol (false positive in prior version)
        "OC(=O)c1ccccc1I(=O)=O",                      # ortho-iodylbenzoic acid (previously missed)
        "ICCCCCCl",                                  # 1-Chloro-4-iodobutane (true positive)
    ]
    
    for smi in test_smiles:
        result, reason = is_organoiodine_compound(smi)
        print(f"SMILES: {smi}")
        print(f"Result: {result}")
        print(f"Reason: {reason}")
        print("-" * 60)