"""
Classifies: CHEBI:33857 aromatic primary alcohol
"""
# aromatic_primary_alcohol.py
"""
Classifies: Aromatic primary alcohol
Definition:
  An aromatic primary alcohol is defined as any primary alcohol in which the alcoholic hydroxy group 
  is attached to a carbon which is itself bonded to an aromatic ring.
  
  This implementation uses a SMARTS substructure search for a benzylic alcohol pattern:
  an aromatic atom (c) directly attached to a CH2 group ([CH2]) that in turn is bonded to an OH group ([OH]).
"""

from rdkit import Chem

def is_aromatic_primary_alcohol(smiles: str):
    """
    Determines if a molecule is an aromatic primary alcohol based on its SMILES string.
    An aromatic primary alcohol has an alcoholic -OH group attached to a primary (CH2) group 
    which is bonded to an aromatic ring.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is classified as an aromatic primary alcohol, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for an aromatic primary alcohol.
    # "c[CH2][OH]" describes:
    #   c    => Aromatic atom (usually carbon in an aromatic ring)
    #   [CH2] => A methylene group (primary carbon: two hydrogens)
    #   [OH]  => An hydroxyl group
    pattern = Chem.MolFromSmarts("c[CH2][OH]")
    
    # Search for the substructure in the molecule.
    matches = mol.GetSubstructMatches(pattern)
    if not matches:
        return False, "No aromatic primary alcohol substructure (c[CH2][OH]) found"

    # If a matching substructure is found, it is classified accordingly.
    return True, "Molecule contains an aromatic primary alcohol group"

# Example usage:
if __name__ == "__main__":
    examples = {
        "4-acetoxybenzyl alcohol": "CC(=O)OC1=CC=C(CO)C=C1",
        "3-pyridinemethanol": "C1=CC(=CN=C1)CO",
        "2-methylbenzyl alcohol": "CC1=C(CO)C=CC=C1",
        "4-cyanobenzyl alcohol": "C1(=CC=C(C=C1)C#N)CO"
    }
    for name, smiles in examples.items():
        result, reason = is_aromatic_primary_alcohol(smiles)
        print(f"{name} ({smiles}): {result}, {reason}")