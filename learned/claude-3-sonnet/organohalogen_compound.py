"""
Classifies: CHEBI:17792 organohalogen compound
"""
"""
Classifies: CHEBI:38114 organohalogen compound
A compound containing at least one carbon-halogen bond (where X is a halogen atom).
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_organohalogen_compound(smiles: str):
    """
    Determines if a molecule is an organohalogen compound based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organohalogen compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Remove salts and counterions
    mol = AllChem.RemoveHeavyAtomSalts(mol)
    
    # Check for carbon-halogen bonds
    has_c_x_bond = False
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() in [9, 17, 35, 53]:
            has_c_x_bond = True
            break
        elif a2.GetAtomicNum() == 6 and a1.GetAtomicNum() in [9, 17, 35, 53]:
            has_c_x_bond = True
            break
    
    if has_c_x_bond:
        return True, "Contains at least one carbon-halogen bond"
    else:
        return False, "No carbon-halogen bonds found"

# Example usage
smiles = "CCCBr"
result, reason = is_organohalogen_compound(smiles)
print(f"Is organohalogen compound? {result} ({reason})")