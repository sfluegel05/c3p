"""
Classifies: CHEBI:17792 organohalogen compound
"""
"""
Classifies: Organohalogen Compound
A compound is considered an organohalogen if it contains at least one carbon–halogen bond.
This implementation uses a SMARTS query to detect bonds between any carbon atom and
one of the halogen atoms: F, Cl, Br or I.
"""

from rdkit import Chem

def is_organohalogen_compound(smiles: str):
    """
    Determines if a molecule is an organohalogen compound based on its SMILES string.
    
    A compound is classified as an organohalogen if it contains at least one carbon–halogen bond.
    That is, at least one halogen atom (F, Cl, Br, or I) must be directly bonded to a carbon atom.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is an organohalogen compound, False otherwise.
        str: Reason for classification.
    """
    # Parse the molecule from the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # To avoid misclassification due to the presence of salts or multiple fragments,
    # use only the largest fragment based on heavy atom count.
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if not frags:
        return False, "No fragments found"
    largest = max(frags, key=lambda m: m.GetNumHeavyAtoms())
    
    # Define a SMARTS pattern that matches a carbon (atomic number 6) directly bonded
    # to any halogen atom (F, Cl, Br, I). 
    # The pattern "[#6]-[F,Cl,Br,I]" will match in either bond direction.
    query = Chem.MolFromSmarts("[#6]-[F,Cl,Br,I]")
    if largest.HasSubstructMatch(query):
        return True, "Molecule contains a carbon–halogen bond"
    
    return False, "No carbon–halogen bond found"

# Example usage for testing:
if __name__ == "__main__":
    test_smiles = [
         "C[C@H]1C[C@@H](C)Br",              # True: carbon–Br bond
         "BrN1C(=O)CCC1=O",                   # False: bromine attached to N (N-bromosuccinimide)
         "OC\\C=C\\Cl",                      # True: carbon–Cl bond via sp2 carbon
         "Oc1cccc(F)c1O",                    # True: aromatic ring with C–F bond
         "Brn1ccc2ccccc12",                  # False: bromine attached to nitrogen (1-bromoindole)
         "Cln1c(=O)n(Cl)c(=O)n(Cl)c1=O",      # False: chlorines on nitrogen (triazinone)
         "OC(=O)c1nc(Cl)ccc1Cl",              # True: at least one Cl is directly bound to C
    ]
    for smi in test_smiles:
        result, reason = is_organohalogen_compound(smi)
        print(f"SMILES: {smi}\n  Result: {result}\n  Reason: {reason}\n")