"""
Classifies: CHEBI:167164 mineral nutrient
"""
"""
Classifies: CHEBI:27027 mineral nutrient
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_mineral_nutrient(smiles: str):
    """
    Determines if a molecule is a mineral nutrient based on its SMILES string.
    A mineral nutrient is typically an inorganic compound containing metal cations
    and inorganic or simple organic anions, essential for human metabolism.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a mineral nutrient, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Extended list of nutrient metal cations
    nutrient_metals = {'Ca', 'Mg', 'K', 'Na', 'Ba', 'Cs', 'Zn', 'Fe', 'Al', 'La', 'Sb', 'Pd', 'Mn', 'Cu', 'Co', 'Mo', 'Se', 'Cr'}
    
    # Check for presence of nutrient metal cations
    has_nutrient_metal = any(atom.GetSymbol() in nutrient_metals for atom in mol.GetAtoms())
    if not has_nutrient_metal:
        return False, "No nutrient metal cation found"

    # List of allowed anions (both inorganic and simple organic)
    allowed_anions = {
        # Inorganic
        'O', 'P', 'S', 'N', 'Cl', 'F', 'Si', 'I', 'Br',
        # Organic (only carboxylates and carbonates)
        'C'
    }
    
    # Check for presence of allowed anions
    has_allowed_anion = any(atom.GetSymbol() in allowed_anions for atom in mol.GetAtoms())
    if not has_allowed_anion:
        return False, "No allowed anion found"

    # Check for organic components
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if len(carbon_atoms) > 0:
        # Only allow carbon in specific patterns (carbonates, carboxylates)
        allowed_carbon_patterns = [
            Chem.MolFromSmarts("[CX3](=[OX1])([OX2])"),  # Carbonates, carboxylates
            Chem.MolFromSmarts("[CX3](=[OX1])([OX1])")   # Carbon dioxide-like
        ]
        
        # Check if all carbon atoms are in allowed patterns
        for carbon in carbon_atoms:
            for pattern in allowed_carbon_patterns:
                if mol.GetAtomWithIdx(carbon.GetIdx()).HasSubstructMatch(pattern):
                    break
            else:
                return False, "Contains carbon atoms not in allowed patterns"

    # Additional checks for complex organic structures
    # Reject molecules with too many rotatable bonds or rings
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable > 2:
        return False, "Too many rotatable bonds for mineral nutrient"
    
    n_rings = rdMolDescriptors.CalcNumRings(mol)
    if n_rings > 1:
        return False, "Too many rings for mineral nutrient"

    # Check molecular weight - mineral nutrients are typically small molecules
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 500:
        return False, "Molecular weight too high for mineral nutrient"

    # Check for ionic character (not strictly required for all cases)
    has_cation = any(atom.GetFormalCharge() > 0 for atom in mol.GetAtoms())
    has_anion = any(atom.GetFormalCharge() < 0 for atom in mol.GetAtoms())
    
    # Allow neutral compounds if they contain nutrient metals and allowed anions
    if not (has_cation or has_anion):
        return True, "Neutral compound containing nutrient metal and allowed anions"

    return True, "Contains nutrient metal cations and allowed anions, typical of mineral nutrients"