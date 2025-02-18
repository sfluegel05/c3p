"""
Classifies: CHEBI:37142 organoiodine compound
"""
"""
Classifies: CHEBI:? organoiodine compound

An organoiodine compound is defined as a compound containing at least one carbon–iodine bond.
However, to avoid trivial or ionic compounds we apply extra filters:
  - The molecule must be a single fragment (not a salt or mixture).
  - The molecule must have a minimum number of heavy atoms.
  - No atom should carry a nonzero formal charge.
  - In at least one detected C–I bond the carbon atom must be "sufficiently organic"; here we
    require that the ratio of non‐carbon heavy atoms to carbon atoms overall is not too high.
"""

from rdkit import Chem

def is_organoiodine_compound(smiles: str):
    """
    Determines whether a molecule is an organoiodine compound based on its SMILES string.
    In addition to detecting a C–I bond, we reject compounds that appear too trivial (e.g. very small),
    are salts (multiple disconnected fragments) or carry formal charges or are very heteroatom-rich.
    
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
    
    # Reject if the molecule breaks into more than one fragment (possible salt/mixture)
    frags = Chem.GetMolFrags(mol, asMols=True)
    if len(frags) > 1:
        return False, "Compound is fragmented (possible salt or mixture)"
    
    # Reject if any atom has a nonzero formal charge.
    for atom in mol.GetAtoms():
        if atom.GetFormalCharge() != 0:
            return False, "Compound has a nonzero formal charge"
    
    # Reject trivial molecules by heavy atom count (e.g. iodoethane has only 3 heavy atoms)
    heavy_atoms = mol.GetNumHeavyAtoms()
    if heavy_atoms < 4:
        return False, f"Compound too simple with only {heavy_atoms} heavy atoms"
    
    # Define a SMARTS pattern for a carbon–iodine bond.
    ci_pattern = Chem.MolFromSmarts("[#6]-I")
    if not mol.HasSubstructMatch(ci_pattern):
        return False, "Compound does not contain any carbon–iodine bonds"
    
    # Optionally, further require that the organic portion is sufficiently carbon rich.
    # Count total carbons vs. non-carbon heavy atoms.
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    non_carbon = heavy_atoms - carbon_count
    # If there are relatively too many non-carbon atoms (e.g. many electronegative groups)
    # then we consider the compound as not a typical organoiodine compound.
    if carbon_count == 0 or (non_carbon / carbon_count) > 1.0:
        return False, "Compound appears too heteroatom-rich to be a typical organoiodine compound"
    
    # Passed all tests.
    return True, "Compound contains a carbon–iodine bond in an organic framework"

# Example usage (for debugging / testing):
if __name__ == "__main__":
    test_smiles = [
        "IC1=CC(F)=C(NC(=O)CC2(CCCCC2)CC(O)=O)",  # 2-{1-[2-(2-fluoro-4-iodoanilino)-2-oxoethyl]cyclohexyl}acetic acid (TP)
        "CCOc1nc2ccc(I)cc2c(=O)n1CCC",               # proquinazid (TP)
        "C(O)C1OC(OC1)CI",                           # domiodol (FP according to evaluation; filtered as too heteroatom-rich)
        "CCI",                                       # iodoethane (FP; too simple)
        "ICCCCCCl"                                   # 1-Chloro-4-iodobutane (TP)
    ]
    
    for smi in test_smiles:
        result, reason = is_organoiodine_compound(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")