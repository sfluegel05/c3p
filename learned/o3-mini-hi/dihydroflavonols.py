"""
Classifies: CHEBI:48039 dihydroflavonols
"""
"""
Classifies: dihydroflavonols 
Definition: Any hydroxyflavanone in which a hydroxy group is present at position 3 of the heterocyclic ring.
That is, a flavanone (benzopyran-4-one scaffold), saturated between C2 and C3, with an –OH group at C3.
"""

from rdkit import Chem

def is_dihydroflavonols(smiles: str):
    """
    Determines if a given SMILES string corresponds to a dihydroflavonol.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is a dihydroflavonol, False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse SMILES string into an RDKit Molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for the dihydroflavonol scaffold.
    # This pattern looks for a benzopyranone system where:
    #   - there is an aromatic ring (ring A) fused to a heterocycle (ring C)
    #   - the heterocycle has an oxygen atom, a carbonyl function (C(=O))
    #     and two saturated (sp3) carbons (CX4) where the first carries an OH (at C3)
    # Note: "c1ccc2c(c1)OC(=O)[CX4](O)[CX4]2" is one possible representation.
    dihydroflavonol_smarts = "c1ccc2c(c1)OC(=O)[CX4](O)[CX4]2"
    pattern = Chem.MolFromSmarts(dihydroflavonol_smarts)
    
    if pattern is None:
        return None, None   # In a case of internal error in pattern
    
    # Check for the dihydroflavonol substructure.
    # This requires the presence of the benzopyran-4-one scaffold with a hydroxyl at the saturated C3 position.
    if mol.HasSubstructMatch(pattern):
        return True, "Found dihydroflavonol scaffold with an OH at the heterocycle's 3-position"
    else:
        return False, "Molecule does not contain the required dihydroflavonol scaffold"

# Example usage (for testing purposes):
if __name__ == "__main__":
    # List of example SMILES from the provided dihydroflavonol examples
    example_smiles = [
        "O1C(C(O)C(=O)C=2C1=CC=3OCOC3C2O",  # 3,5-Dihydroxy-6,7-methylenedioxyflavanone (note: missing closing parentheses may be present in an example)
        "O[C@@H]1[C@H](Oc2cc(O)ccc2C1=O)c1ccc(O)cc1",  # garbanzol
        "OC1C(Oc2cc(O)cc(O)c2C1=O)c1ccc(O)c(O)c1",  # (+)-taxifolin (one of the examples)
    ]
    
    for s in example_smiles:
        result, reason = is_dihydroflavonols(s)
        print(f"SMILES: {s}\nClassified as dihydroflavonol? {result}\nReason: {reason}\n")