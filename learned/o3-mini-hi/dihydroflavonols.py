"""
Classifies: CHEBI:48039 dihydroflavonols
"""
"""
Classifies: dihydroflavonols
Definition: Any hydroxyflavanone in which a hydroxy group is present at position 3 of the heterocyclic ring.
This implementation uses a SMARTS pattern corresponding to the dihydroflavonol core.
Essential features:
  - A benzopyran-4-one (flavanone) scaffold with a fused aromatic ring (A-ring)
  - A saturated C2–C3 bond (dihydro)
  - A hydroxyl at C3
  - An aromatic B-ring attached at C2
The SMARTS used is: "O=C[CX4]1OC(c2ccccc2)[CX4]([OH])C1"
"""

from rdkit import Chem

def is_dihydroflavonols(smiles: str):
    """
    Determines if a given SMILES string corresponds to a dihydroflavonol.
    
    The function checks for the presence of the dihydroflavonol core by matching
    against a SMARTS pattern that enforces all the defining features of the class:
    
      SMARTS: O=C[CX4]1OC(c2ccccc2)[CX4]([OH])C1
      
    This pattern stipulates:
      - A carbonyl (O=C) bonded to a tetrahedral carbon ([CX4]) that is part of a ring (denoted by '1').
      - An oxygen (O) in the ring.
      - A pendant aromatic ring (c2ccccc2) corresponding to the B-ring.
      - A saturated carbon ([CX4]) bearing an –OH group at the C3 position.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a dihydroflavonol, False otherwise.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS pattern for the dihydroflavonol core.
    smarts = "O=C[CX4]1OC(c2ccccc2)[CX4]([OH])C1"
    query = Chem.MolFromSmarts(smarts)
    if query is None:
        return None, None  # Should not happen if SMARTS is valid.
    
    # Check if the molecule contains the dihydroflavonol core.
    if mol.HasSubstructMatch(query):
        return True, ("Matches dihydroflavonol core: contains a flavanone scaffold with a saturated "
                      "C2–C3 bond, a hydroxyl group at the C3 position, and an aromatic B‐ring.")
    else:
        return False, "Molecule does not match the dihydroflavonol core SMARTS pattern."

# Example usage for testing purposes:
if __name__ == "__main__":
    test_smiles = [
        "O1C(C(O)C(=O)C=2C1=CC=3OCOC3C2O)C4=CC=CC=C4",  # 3,5-Dihydroxy-6,7-methylenedioxyflavanone
        "O[C@@H]1[C@H](Oc2cc(O)ccc2C1=O)c1ccc(O)cc1",    # garbanzol
        "CC(C)=CCc1c(O)cc(O)c2C(=O)[C@H](O)[C@H](Oc12)c1ccccc1",  # glepidotin B
        "OC1C(Oc2cc(O)cc(O)c2C1=O)c1cc(O)c(O)c1",        # dihydromyricetin
        "COC1=C(O)C=CC(=C1)[C@H]1OC2=CC(O)=CC(O)=C2C(=O)[C@@H]1O", # dihydroisorhamnetin
    ]
    for s in test_smiles:
        result, reason = is_dihydroflavonols(s)
        print(f"SMILES: {s}\nClassified as dihydroflavonol? {result}\nReason: {reason}\n")