"""
Classifies: CHEBI:26658 sesquiterpenoid
"""
"""
Classifies: Sesquiterpenoid (any terpenoid derived from a sesquiterpene with a C15 skeleton that may be rearranged or modified)
"""
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import rdMolDescriptors

def is_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a sesquiterpenoid based on its SMILES string.
    Sesquiterpenoids are derived from sesquiterpenes (C15 skeleton) but may have rearrangements or modifications.
    This heuristic extracts the Murcko scaffold (the core scaffold) and counts the number of carbon atoms.
    If the core has roughly 15 carbons (allowing a window of 13 to 16 to tolerate minor modifications),
    the molecule is considered a sesquiterpenoid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is classified as a sesquiterpenoid, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string to produce an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Extract the Murcko scaffold (the molecular framework) of the molecule.
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    if scaffold is None:
        return False, "Unable to extract scaffold from the molecule"
        
    # Count the number of carbon atoms present in the scaffold.
    scaffold_c_count = sum(1 for atom in scaffold.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # A sesquiterpene skeleton has 15 carbons. Allow slight variation for rearrangements (e.g., 13–16).
    if 13 <= scaffold_c_count <= 16:
        return True, f"Scaffold contains {scaffold_c_count} carbons, consistent with a sesquiterpene-derived skeleton."
    else:
        return False, f"Scaffold contains {scaffold_c_count} carbons, not consistent with the expected ~15 carbons of a sesquiterpene-derived structure."

# Example usage (for testing purposes):
if __name__ == '__main__':
    test_smiles = "CC(=C)[C@@H]1CCC(CO)=C2CC(=O)C(C)=C2C1"  # Example (indicanone, may or may not qualify depending on the scaffold)
    result, reason = is_sesquiterpenoid(test_smiles)
    print(result, reason)