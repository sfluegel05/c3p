"""
Classifies: CHEBI:27933 beta-lactam antibiotic
"""
"""
Classifies: Beta-lactam antibiotics
Definition: An organonitrogen heterocyclic antibiotic that contains a beta-lactam ring.
A beta-lactam ring is defined as a four-membered cyclic amide containing one nitrogen,
three carbons, and one of the carbons bears a carbonyl group (C=O).
In this implementation we search for a substructure matching the SMARTS pattern:
    N1C(=O)[C;R][C;R]1
which looks for a 4-membered ring (all three carbon atoms must be in a ring) with one nitrogen and a carbonyl.
Note: This is only a structural filter and may (or may not) include edge cases.
"""

from rdkit import Chem

def is_beta_lactam_antibiotic(smiles: str):
    """
    Determines if a molecule is a beta-lactam antibiotic based on its SMILES string.
    We require that the molecule contains at least one beta-lactam ring:
      a four-membered cyclic amide having one nitrogen, three carbons,
      and with one of the carbon atoms bearing a double bond to oxygen.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule appears to contain a beta-lactam ring, False otherwise.
        str: Explanation for the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # An issue can occur if the bond orders are not fully explicit.
    # We try to kekulize the molecule to normalize bond orders.
    try:
        Chem.Kekulize(mol, clearAromaticFlags=True)
    except Exception:
        # If kekulization fails, we continue anyway.
        pass

    # Define a SMARTS pattern for a beta-lactam ring.
    # This pattern requires:
    #  - A ring numbered 1
    #  - A nitrogen (N1)
    #  - A carbon carrying a double bond to an oxygen: C(=O)
    #  - Two additional ring carbons (specified as [C;R] to be in a ring)
    beta_lactam_smarts = Chem.MolFromSmarts("N1C(=O)[C;R][C;R]1")
    if beta_lactam_smarts is None:
        return False, "Error creating SMARTS pattern"
    
    # Look for a matching substructure.
    if mol.HasSubstructMatch(beta_lactam_smarts):
        return True, ("Molecule contains a beta-lactam ring "
                      "(a 4-membered cyclic amide with 1 N, 3 C, and one C=O), "
                      "which is typical for beta-lactam antibiotics.")
    else:
        return False, "No beta-lactam ring found in the molecule"

# (Optional) Simple testing – these examples are for illustration.
if __name__ == "__main__":
    # List a couple of example SMILES strings.
    smiles_examples = [
        "[H][C@]12SCC(C)=C(N1C(=O)[C@H]2N)C(O)=O",    # 7beta-aminodeacetoxycephalosporanic acid (true positive)
        "O=C1CCN1",                                   # azetidin-2-one (true positive)
        "C1=CC=CC=C1"                                  # benzene (should be false)
    ]
    for smi in smiles_examples:
        res, reason = is_beta_lactam_antibiotic(smi)
        print(f"SMILES: {smi}\nResult: {res}\nReason: {reason}\n")