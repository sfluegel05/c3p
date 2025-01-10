"""
Classifies: CHEBI:37142 organoiodine compound
"""
"""
Classifies: CHEBI:24863 organoiodine compound
Definition: A compound containing at least one carbon-iodine bond.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_organoiodine_compound(smiles: str):
    """
    Determines if a molecule is an organoiodine compound based on its SMILES string.
    An organoiodine compound must contain at least one carbon-iodine bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organoiodine compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find all iodine atoms
    iodine_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 53]
    
    if not iodine_atoms:
        return False, "No iodine atoms found"

    # Check for valid carbon-iodine covalent bonds
    for iodine in iodine_atoms:
        # Skip if iodine has formal charge
        if iodine.GetFormalCharge() != 0:
            continue
            
        for bond in iodine.GetBonds():
            # Get the atom at the other end of the bond
            other_atom = bond.GetOtherAtom(iodine)
            
            # Check if it's a carbon atom
            if other_atom.GetAtomicNum() != 6:
                continue
                
            # Check if carbon has reasonable valence and no charge
            if other_atom.GetFormalCharge() != 0:
                continue
                
            # Verify it's a single covalent bond
            if bond.GetBondType() != Chem.BondType.SINGLE:
                continue
                
            # Additional checks for reasonable geometry
            if other_atom.GetDegree() > 4:
                continue
                
            return True, "Contains at least one carbon-iodine covalent bond"
    
    return False, "No valid carbon-iodine covalent bonds found"