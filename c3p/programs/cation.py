"""
Classifies: CHEBI:36916 cation
"""
"""
Classifies chemical entities of the class cation:
A monoatomic or polyatomic species having one or more elementary charges of the proton.
This improved classifier not only checks the overall net formal charge,
but also inspects for permanent cationic functional groups (e.g., ammonium groups)
that may be present in zwitterionic species.
"""

from rdkit import Chem

def is_cation(smiles: str):
    """
    Determines if a molecule is a cation based on its SMILES string.
    
    A molecule is classified as a cation if:
      - It has a net formal charge greater than zero,
    OR
      - It has an overall net zero charge but contains a permanent cationic functional group
        (such as an ammonium group) indicating the presence of a protonated moiety.
        
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a cation, False otherwise.
        str: Explanation of the classification.
    """
    # Parse the SMILES string to create a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Calculate net formal charge by summing formal charges of all atoms
    net_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    
    # If overall net charge is > 0 then it is a cation.
    if net_charge > 0:
        return True, f"Molecule has a net positive charge of {net_charge}"
    
    # When net charge is zero, the molecule might be zwitterionic.
    # Look for permanent cationic functional groups (e.g. ammonium groups).
    # We define a list of SMARTS patterns for common cationic moieties.
    cation_smarts = [
        "[NH3+]",       # protonated ammonium with three hydrogens
        "[N+](C)(C)C",  # typical trimethylammonium group
        "[NH2+]",       # protonated amine with two hydrogens
        "[N+;H]"        # general pattern: nitrogen with at least one hydrogen and a positive charge
    ]
    
    for smarts in cation_smarts:
        patt = Chem.MolFromSmarts(smarts)
        if patt is not None and mol.HasSubstructMatch(patt):
            return True, "Molecule contains a permanent cationic functional group despite net zero charge"
    
    # If net charge is negative or zero and no cationic group was found, then it is not a cation.
    if net_charge < 0:
        return False, f"Molecule has a net negative charge of {net_charge} (not a cation)"
    else:
        return False, f"Molecule has a net charge of {net_charge} and no identified permanent cationic group"
        
# Example usage (uncomment the following lines to test):
# smiles_examples = [
#     "C[N+](C)(C)C",  # trimethylammonium: should be a cation (net +1)
#     "OCC[NH+](C)C",  # protonated amine: cation (net +1)
#     "P(OCC[N+](C)(C)C)([O-])=O",  # zwitterionic phosphocholine: net 0 but contains [N+](C)(C)C -> cation
#     "[O-]C(=O)C",    # acetate anion (net -1): not a cation
#     "CCO"            # ethanol (net 0, no cationic group): not a cation
# ]
#
# for smi in smiles_examples:
#     result, reason = is_cation(smi)
#     print(f"SMILES: {smi}\nResult: {result} | Reason: {reason}\n")