"""
Classifies: CHEBI:37142 organoiodine compound
"""
"""
Classifies: CHEBI:? organoiodine compound

An organoiodine compound is defined as a compound containing at least one carbon–iodine bond.
Our improved classifier no longer restricts iodine to a monovalent bonding environment so that
hypervalent (e.g. iodyl) iodine compounds are also captured. (Optionally one could filter out very
simple compounds by, e.g., checking the total heavy-atom count.)
"""

from rdkit import Chem

def is_organoiodine_compound(smiles: str):
    """
    Determines if a molecule is an organoiodine compound based on its SMILES string.
    The definition is "compound contains at least one carbon–iodine bond" regardless of whether
    the iodine atom carries additional substituents (for example I(=O)=O groups).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule contains a carbon–iodine bond, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string to get the molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Optionally: filter out very trivial (low molecular weight) compounds.
    # Uncomment the following lines if you wish to require, for example, at least 4 heavy atoms.
    # heavy_atoms = mol.GetNumHeavyAtoms()
    # if heavy_atoms < 4:
    #     return False, f"Compound too simple ({heavy_atoms} heavy atoms) to be an organoiodine compound."
    
    # Define a SMARTS pattern that matches any bond between a carbon (atomic num 6) and iodine (atomic num 53).
    ci_pattern = Chem.MolFromSmarts("[#6]-I")
    
    # Check if the molecule contains at least one carbon–iodine bond.
    if mol.HasSubstructMatch(ci_pattern):
        return True, "Compound contains at least one carbon–iodine bond."
    else:
        return False, "Compound does not contain any carbon–iodine bonds."

# Example usage:
if __name__ == "__main__":
    test_smiles = [
        "IC1=CC(F)=C(NC(=O)CC2(CCCCC2)CC(O)=O)",  # 2-{1-[2-(2-fluoro-4-iodoanilino)-2-oxoethyl]cyclohexyl}acetic acid (TP)
        "CCOc1nc2ccc(I)cc2c(=O)n1CCC",              # proquinazid (TP)
        "C(O)C1OC(OC1)CI",                          # domiodol (FP according to evaluation that considers this too trivial)
        "CCI"                                       # iodoethane (FP: simple alkyl iodide)
    ]
    for smi in test_smiles:
        result, reason = is_organoiodine_compound(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")