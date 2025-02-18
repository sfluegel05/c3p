"""
Classifies: CHEBI:17792 organohalogen compound
"""
"""
Classifies: organohalogen compound (CHEBI: orgHalogen)
"""
from rdkit import Chem

def is_organohalogen_compound(smiles: str):
    """
    Determines if a molecule is an organohalogen compound by checking for at least one carbon-halogen bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains a carbon-halogen bond, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Halogens: F, Cl, Br, I, At
    halogens = {9, 17, 35, 53, 85}
    halogen_symbols = {9: 'F', 17: 'Cl', 35: 'Br', 53: 'I', 85: 'At'}
    
    # Iterate through all bonds to check for C-X bonds
    for bond in mol.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        atoms = [atom1, atom2]
        
        # Check if the bond connects a carbon and a halogen
        has_carbon = any(a.GetAtomicNum() == 6 for a in atoms)
        has_halogen = any(a.GetAtomicNum() in halogens for a in atoms)
        
        if has_carbon and has_halogen:
            # Identify the halogen atom
            halogen = next(a for a in atoms if a.GetAtomicNum() in halogens)
            symbol = halogen_symbols[halogen.GetAtomicNum()]
            return True, f"Contains carbon-{symbol} bond"
    
    return False, "No carbon-halogen bonds found"