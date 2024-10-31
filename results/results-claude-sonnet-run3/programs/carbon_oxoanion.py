from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors

def is_carbon_oxoanion(smiles: str):
    """
    Determines if a molecule is a carbon oxoanion (negative ion containing only C and O).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carbon oxoanion, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Check formal charge
    total_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    if total_charge >= 0:
        return False, "Not an anion (total charge >= 0)"

    # Get all atoms
    atoms = mol.GetAtoms()
    
    # Check if molecule contains only C and O atoms
    allowed_atoms = {'C', 'O'}
    atom_symbols = {atom.GetSymbol() for atom in atoms}
    if not atom_symbols.issubset(allowed_atoms):
        disallowed = atom_symbols - allowed_atoms
        return False, f"Contains atoms other than C and O: {disallowed}"

    # Count C and O atoms
    c_count = len([atom for atom in atoms if atom.GetSymbol() == 'C'])
    o_count = len([atom for atom in atoms if atom.GetSymbol() == 'O'])
    
    # Count negative charges
    neg_charge_count = len([atom for atom in atoms if atom.GetFormalCharge() < 0])
    
    # Analyze the structure
    if c_count == 0:
        return False, "No carbon atoms present"
        
    if o_count == 0:
        return False, "No oxygen atoms present"
        
    if neg_charge_count == 0:
        return False, "No negative charges present"

    return True, f"Carbon oxoanion with formula C{c_count}O{o_count}({total_charge}-)"
# Pr=1.0
# Recall=0.8571428571428571