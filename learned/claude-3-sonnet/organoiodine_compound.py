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

    # Check for radioactive isotopes or charged species
    for iodine in iodine_atoms:
        if iodine.GetIsotope() > 0:  # Skip radioactive isotopes
            continue
        if iodine.GetFormalCharge() != 0:  # Skip charged iodine
            continue
            
        # Check bonds to iodine
        for bond in iodine.GetBonds():
            # Get the atom at the other end of the bond
            other_atom = bond.GetOtherAtom(iodine)
            
            # Must be a carbon atom
            if other_atom.GetAtomicNum() != 6:
                continue
                
            # Carbon must be neutral
            if other_atom.GetFormalCharge() != 0:
                continue
                
            # Must be a single covalent bond
            if bond.GetBondType() != Chem.BondType.SINGLE:
                continue
                
            # Carbon must have valid valence
            if other_atom.GetTotalValence() > 4:
                continue
                
            # Carbon should not be part of charged/ionic groups
            valid_carbon = True
            for neighbor in other_atom.GetNeighbors():
                if neighbor.GetFormalCharge() != 0:
                    valid_carbon = False
                    break
                    
            if not valid_carbon:
                continue
                
            # Check if this is part of a salt/ionic compound
            fragment = Chem.rdmolops.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
            if len(fragment) > 1:
                # If molecule has multiple fragments, ensure the C-I bond is in an organic fragment
                for frag in fragment:
                    if other_atom.GetIdx() in [a.GetIdx() for a in frag.GetAtoms()]:
                        # Count C, H, O, N in fragment to verify it's organic
                        c_count = sum(1 for a in frag.GetAtoms() if a.GetAtomicNum() == 6)
                        if c_count < 1:
                            continue
            
            return True, "Contains at least one carbon-iodine covalent bond"
    
    return False, "No valid carbon-iodine covalent bonds found"