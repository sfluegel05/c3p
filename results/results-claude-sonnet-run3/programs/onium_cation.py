from rdkit import Chem
from rdkit.Chem import AllChem

def is_onium_cation(smiles: str):
    """
    Determines if a molecule is an onium cation (mononuclear cation derived by addition of a hydron 
    to a mononuclear parent hydride of pnictogen, chalcogen or halogen families)

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an onium cation, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define valid central atoms for onium cations
    pnictogens = ['N', 'P', 'As', 'Sb', 'Bi']  # Group 15
    chalcogens = ['O', 'S', 'Se', 'Te']         # Group 16
    halogens = ['F', 'Cl', 'Br', 'I']           # Group 17
    valid_central_atoms = set(pnictogens + chalcogens + halogens)

    # Count atoms by type
    atom_counts = {}
    charged_atoms = []
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        atom_counts[symbol] = atom_counts.get(symbol, 0) + 1
        if atom.GetFormalCharge() > 0:
            charged_atoms.append(atom)

    # Should have exactly one positively charged atom
    if len(charged_atoms) != 1:
        return False, "Must have exactly one positively charged atom"
    
    central_atom = charged_atoms[0]
    if central_atom.GetFormalCharge() != 1:
        return False, "Charge must be +1"

    # Central atom must be from valid groups
    symbol = central_atom.GetSymbol()
    if symbol not in valid_central_atoms:
        return False, f"{symbol} is not a valid central atom for onium cations"

    # All other atoms must be hydrogens
    for atom_symbol in atom_counts:
        if atom_symbol != 'H' and atom_symbol != symbol:
            return False, f"Contains non-hydrogen atom: {atom_symbol}"

    # Get total hydrogens including explicit and implicit
    h_count = len([atom for atom in central_atom.GetNeighbors() if atom.GetSymbol() == 'H'])
    h_count += central_atom.GetNumImplicitHs()

    # Check number of hydrogens based on group
    if symbol in pnictogens and h_count != 4:
        return False, f"Pnictogen onium cation should have 4 hydrogens"
    elif symbol in chalcogens and h_count != 3:
        return False, f"Chalcogen onium cation should have 3 hydrogens"
    elif symbol in halogens and h_count != 2:
        return False, f"Halogen onium cation should have 2 hydrogens"

    # Classify based on central atom
    onium_names = {
        'N': 'ammonium',
        'P': 'phosphonium', 
        'As': 'arsonium',
        'Sb': 'stibonium',
        'Bi': 'bismuthonium',
        'O': 'oxonium',
        'S': 'sulfonium',
        'Se': 'selenonium',
        'Te': 'telluronium',
        'F': 'fluoronium',
        'Cl': 'chloronium', 
        'Br': 'bromonium',
        'I': 'iodonium'
    }

    return True, f"Valid {onium_names[symbol]} cation"
# Pr=None
# Recall=0.0