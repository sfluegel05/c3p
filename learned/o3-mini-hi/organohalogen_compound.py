"""
Classifies: CHEBI:17792 organohalogen compound
"""
"""
Classifies: Organohalogen Compound
A compound is considered an organohalogen if it contains at least one carbon–halogen bond.
This implementation parses the molecule (using only its largest fragment to avoid salts), then
iterates over each halogen atom and checks if any bonded neighbor is carbon.
"""

from rdkit import Chem

def is_organohalogen_compound(smiles: str):
    """
    Determines if a molecule is an organohalogen compound based on its SMILES string.
    
    A compound is classified as an organohalogen if it contains at least one 
    carbon–halogen bond (halogen = F, Cl, Br or I). That is, at least one halogen atom
    must be directly bonded to a carbon atom.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is an organohalogen compound, False otherwise.
        str: Reason for classification.
    """
    # Parse molecule from SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Use only the largest fragment to avoid counting counterions
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if not frags:
        return False, "No fragments found"
    largest = max(frags, key=lambda m: m.GetNumHeavyAtoms())
    
    # Iterate over all atoms in the largest fragment.
    for atom in largest.GetAtoms():
        # If atom is a halogen
        if atom.GetSymbol() in {"F", "Cl", "Br", "I"}:
            # Check if any neighbor is a carbon atom (atomic number 6)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6:
                    return True, "Molecule contains a carbon–halogen bond"
    
    # If no halogen directly bonded to a carbon was found:
    return False, "No carbon–halogen bond found"

# Example usage (for testing)
if __name__ == "__main__":
    test_smiles = [
         "C[C@H]1C[C@@H](C)Br",             # Direct C–Br bond
         "BrN1C(=O)CCC1=O",                  # Bromine attached only to N (not a C–X bond)
         "OC\\C=C\\Cl",                     # Direct C–Cl bond via an sp2 carbon
         "Oc1cccc(F)c1O",                   # Aromatic ring with C–F bond
         "Brn1ccc2ccccc12",                 # 1-bromoindole; Br attached to N => NOT organohalogen
         "Cln1c(=O)n(Cl)c(=O)n(Cl)c1=O",     # Trichloro triazinone; chlorines on N => NOT organohalogen
         "OC(=O)c1nc(Cl)ccc1Cl",             # clopyralid; at least one Cl is directly bound to C
    ]
    for smi in test_smiles:
        result, reason = is_organohalogen_compound(smi)
        print(f"SMILES: {smi}\n  Result: {result}\n  Reason: {reason}\n")