"""
Classifies: CHEBI:28863 flavanones
"""
#!/usr/bin/env python3
"""
Classifies: flavanones – Members of the class of flavans with a 
3,4-dihydro-2-aryl-2H-1-benzopyran-4-one skeleton and its substituted derivatives.
This improved pattern detects the bicyclic chroman-4-one core and enforces that the C-2 (saturated) 
carbon bears an aromatic substituent.
"""
from rdkit import Chem

def is_flavanones(smiles: str):
    """
    Determines if a molecule is a flavanone based on its SMILES string.
    
    Flavanones are defined by having the 3,4-dihydro-2-aryl-2H-1-benzopyran-4-one core.
    This implementation uses a SMARTS pattern that looks for a fused bicyclic system in which:
      - The chroman-4-one (benzopyran-4-one) skeleton is present (an oxygen-containing ring fused 
        to a benzene ring with a ketone function) 
      - The saturated (C-2) position bears an aromatic substituent ([a]), which is our proxy for the “2-aryl” part.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is annotated as a flavanone, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a refined SMARTS pattern representing the flavanone core.
    # Explanation of the pattern "c1ccc2OC(C([a])C2=O)c1":
    #   - "c1ccc2"        : An aromatic ring (ring A) that is fused to ring 2.
    #   - "OC("           : The fused ring 2 starts with an oxygen.
    #   - "C([a])"        : The next (saturated) carbon (position 2) bears a substituent that
    #                       is forced to be an aromatic atom ([a]) indicating the aryl group.
    #   - "C2=O)"        : Closes the 6-membered dihydropyran ring by a carbonyl (the 4-one).
    #   - "c1"            : Completes the fusion back into the A ring.
    core_smarts = "c1ccc2OC(C([a])C2=O)c1"
    core_query = Chem.MolFromSmarts(core_smarts)
    if core_query is None:
        return False, "Error in SMARTS pattern for flavanone core"
    
    # Search for the flavanone core in the molecule (ignoring chirality)
    matches = mol.GetSubstructMatches(core_query, useChirality=False)
    if matches:
        return True, "Molecule contains the flavanone core (3,4-dihydro-2-aryl-2H-1-benzopyran-4-one)"
    else:
        return False, "Flavanone core structure (3,4-dihydro-2-aryl-2H-1-benzopyran-4-one) not found"

# Example usage (for testing purposes)
if __name__ == "__main__":
    # Some example SMILES strings for flavanones (true positives) and a few known negatives.
    examples = {
        "eriodictyol": "Oc1cc(O)c2C(=O)C[C@H](Oc2c1)c1ccc(O)c(O)c1",
        "sakuranetin": "COc1cc(O)c2C(=O)C[C@H](Oc2c1)c1ccc(O)cc1",
        "hesperetin": "COc1ccc(cc1O)[C@@H]1CC(=O)c2c(O)cc(O)cc2O1",
        "(R)-naringenin": "Oc1ccc(cc1)[C@H]1CC(=O)c2c(O)cc(O)cc2O1",
        # A false positive example from previous trials (sulfonated structure)
        "false_positive": "S(OC=1C(O)=CC(C2OC3=C(C(=O)C2O)C(O)=C(O)C(O)=C3)=CC1O)(O)(=O)=O",
        # A false negative example (diosmin) – should be recognized as a flavanone derivative
        "diosmin": "COc1ccc(cc1O)-c1cc(=O)c2c(O)cc(O)c2o1"
    }
    
    for name, smi in examples.items():
        result, reason = is_flavanones(smi)
        print(f"{name}: {result} ({reason})")