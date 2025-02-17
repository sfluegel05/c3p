"""
Classifies: CHEBI:25409 monoterpenoid
"""
"""
Classifies: Monoterpenoid (any terpenoid derived from a monoterpene, C10 skeleton possibly rearranged or modified)
Heuristic:
  - First, parse the molecule from SMILES.
  - Derive the Bemis–Murcko scaffold (which should contain the underlying carbon skeleton).
  - Count the number of carbons in the scaffold.
  - If the number of carbons is roughly in the range expected for a monoterpene-derived skeleton (here chosen as 7 to 13),
    then we classify the molecule as a monoterpenoid. Otherwise, it is not.
Note: This is a heuristic classification.
"""
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

def is_monoterpenoid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid based on its SMILES string.
    A monoterpenoid is defined as any terpenoid derived from a monoterpene. The criteria used here are based
    on extracting the Bemis–Murcko scaffold and checking whether its carbon count falls in a typical range.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        (bool, str): Tuple where the first element indicates if the molecule is classified as a monoterpenoid,
                     and the second element is the reason for the classification decision.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    try:
        # Obtain the Bemis–Murcko scaffold of the molecule.
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    except Exception as e:
        return False, f"Error generating Murcko scaffold: {e}"
    
    if scaffold is None:
        return False, "Could not generate a scaffold."

    # Count the number of carbon atoms in the scaffold.
    scaffold_carbon_count = sum(1 for atom in scaffold.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Heuristic check: For monoterpenoids, one would expect a core scaffold close to a C10 structure.
    # Allowing for slight variations (loss of a methyl or addition of rearrangements), we choose an inclusive range.
    if scaffold_carbon_count < 7:
        return False, f"Scaffold carbon count too low ({scaffold_carbon_count}). Expected at least 7 carbons."
    elif scaffold_carbon_count > 13:
        return False, f"Scaffold carbon count too high ({scaffold_carbon_count}). Expected 13 or fewer carbons."
    
    return True, f"Scaffold carbon count is {scaffold_carbon_count}, consistent with a monoterpenoid core."