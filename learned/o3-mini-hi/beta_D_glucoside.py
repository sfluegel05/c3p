"""
Classifies: CHEBI:22798 beta-D-glucoside
"""
"""
Classifies: beta-D-glucoside
Definition: Any D-glucoside in which the anomeric centre has beta-configuration.
This routine attempts to detect a beta-D-glucoside moiety by searching for a pyranose ring
substructure attached via an oxygen to an aglycone, where the anomeric centre (the first chiral
carbon of the ring) has beta stereochemistry (encoded as @@ in the SMILES).
"""

from rdkit import Chem

def is_beta_D_glucoside(smiles: str):
    """
    Determines if a molecule is a beta-D-glucoside based on its SMILES string.
    We classify a molecule as a beta-D-glucoside if it contains a pyranose (six-membered sugar) ring
    attached through an oxygen that has beta stereochemistry at the anomeric centre.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a beta-D-glucoside, False otherwise.
        str : Reason for the classification.
    """
    
    # Parse SMILES to an RDKit Mol object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS for beta-D-glucoside.
    # Explanation:
    #   - [!H0] forces that the oxygen is attached to a heavy atom (the aglycone).
    #   - O[C@@H]1 ... indicates that the anomeric carbon (C1) is drawn with the @@ stereochemistry,
    #     which in our examples denotes beta configuration.
    #   - The rest of the pattern [O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O] attempts to capture the
    #     pyranose ring of D-glucose.
    glucoside_smarts = "[!H0]-O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O"
    glucoside_pattern = Chem.MolFromSmarts(glucoside_smarts)
    if glucoside_pattern is None:
        return False, "Error in SMARTS pattern definition"
    
    # Search for the beta-D-glucoside substructure within the molecule.
    if mol.HasSubstructMatch(glucoside_pattern):
        return True, "Molecule contains a beta-D-glucoside moiety with beta-anomeric configuration."
    else:
        return False, "No beta-D-glucoside substructure (with beta anomeric centre) found."

# Example test (uncomment to run)
# smiles_example = "CCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)[C@H](O)\\C=C\\CCCCCCCCCCCCC"
# result, reason = is_beta_D_glucoside(smiles_example)
# print(result, reason)