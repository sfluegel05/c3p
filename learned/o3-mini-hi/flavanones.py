"""
Classifies: CHEBI:28863 flavanones
"""
#!/usr/bin/env python3
"""
Classifies: flavanones – Members of the class of flavans with a 
3,4-dihydro-2-aryl-2H-1-benzopyran-4-one skeleton and its substituted derivatives.
This pattern is detected by searching for a bicyclic benzopyran-4-one core 
with an extra aryl substituent at the C-2 position.
The SMARTS pattern used is: "c1ccc2C(=O)CC([$([c])])Oc2c1"
which represents the 2-phenylchroman-4-one scaffold.
"""
from rdkit import Chem

def is_flavanones(smiles: str):
    """
    Determines if a molecule is a flavanone based on its SMILES string.
    
    Flavanones are defined by the presence of the 3,4-dihydro-2-aryl-2H-1-benzopyran-4-one scaffold.
    In this approach we search for a benzopyran-4-one fused system where:
      - A dihydropyran ring (containing one oxygen and a ketone at the 4 position)
      - Is fused to an aromatic ring (the A ring)
      - And has an extra aryl substituent attached at the C-2 position (verified by forcing the substituent 
        to start with an aromatic atom via the SMARTS constraint [$([c])]).
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a flavanone, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Define a SMARTS pattern representing the 2-phenylchroman-4-one core.
    # Explanation:
    #  "c1ccc2"      : an aromatic ring (ring 1) with an open connection to ring 2.
    #  "C(=O)"       : a ketone function in ring 2 (C-4 position).
    #  "CC([$([c])])": two saturated carbons in ring 2; the second bears a substituent;
    #                   that substituent is forced to begin with an aromatic atom [$([c])],
    #                   which is our probe for the extra (2-aryl) group.
    #  "Oc2c1"       : closing the fused ring by an oxygen and connecting back to the aromatic ring.
    core_smarts = "c1ccc2C(=O)CC([$([c])])Oc2c1"
    core_query = Chem.MolFromSmarts(core_smarts)
    if core_query is None:
        return False, "Error in SMARTS pattern for flavanone core"
    
    # Search for the flavanone core in the molecule (ignore chirality).
    core_matches = mol.GetSubstructMatches(core_query, useChirality=False)
    if core_matches:
        return True, "Molecule contains the flavanone core (3,4-dihydro-2-aryl-2H-1-benzopyran-4-one)"
    else:
        return False, "Flavanone core structure (3,4-dihydro-2-aryl-2H-1-benzopyran-4-one) not found"

# Example usage (for testing purposes)
if __name__ == "__main__":
    examples = {
        "eriodictyol": "Oc1cc(O)c2C(=O)C[C@H](Oc2c1)c1ccc(O)c(O)c1",
        "sakuranetin": "COc1cc(O)c2C(=O)C[C@H](Oc2c1)c1ccc(O)cc1",
        "hesperetin": "COc1ccc(cc1O)[C@@H]1CC(=O)c2c(O)cc(O)cc2O1",
        "(R)-naringenin": "Oc1ccc(cc1)[C@H]1CC(=O)c2c(O)cc(O)cc2O1"
    }
    for name, smi in examples.items():
        result, reason = is_flavanones(smi)
        print(f"{name}: {result} ({reason})")