"""
Classifies: CHEBI:46662 mineral
"""
"""
Classifies: CHEBI:33686 mineral
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_mineral(smiles: str):
    """
    Determines if a molecule is a mineral based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a mineral, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get all atoms
    atoms = mol.GetAtoms()
    if len(atoms) == 0:
        return False, "No atoms found"
    
    # Define common elements in minerals
    mineral_metals = {11: 'Na', 12: 'Mg', 13: 'Al', 19: 'K', 20: 'Ca', 
                     26: 'Fe', 28: 'Ni', 29: 'Cu', 30: 'Zn', 56: 'Ba', 
                     55: 'Cs', 51: 'Sb', 33: 'As'}
    
    mineral_nonmetals = {8: 'O', 9: 'F', 15: 'P', 16: 'S', 17: 'Cl', 
                        14: 'Si', 5: 'B', 7: 'N'}
    
    # Track atom types and properties
    metal_ions = set()
    nonmetal_ions = set()
    has_metal = False
    c_count = 0
    c_c_bonds = 0
    charged_atoms = []
    aromatic_atoms = 0
    
    # Count bonds between carbons
    for bond in mol.GetBonds():
        if (bond.GetBeginAtom().GetAtomicNum() == 6 and 
            bond.GetEndAtom().GetAtomicNum() == 6):
            c_c_bonds += 1
    
    # Analyze atoms
    for atom in atoms:
        atomic_num = atom.GetAtomicNum()
        formal_charge = atom.GetFormalCharge()
        
        if atom.GetIsAromatic():
            aromatic_atoms += 1
            
        if atomic_num in mineral_metals:
            has_metal = True
            metal_ions.add(mineral_metals[atomic_num])
            if formal_charge != 0:
                charged_atoms.append(formal_charge)
                
        elif atomic_num in mineral_nonmetals:
            nonmetal_ions.add(mineral_nonmetals[atomic_num])
            if formal_charge != 0:
                charged_atoms.append(formal_charge)
                
        elif atomic_num == 6:  # Carbon
            c_count += 1
    
    # Immediate exclusion rules
    if aromatic_atoms > 0:
        return False, "Contains aromatic rings - not a mineral"
        
    if c_c_bonds > 2 and not any(x in ["Mg", "Ca", "Ba"] for x in metal_ions):
        return False, "Too many C-C bonds for a mineral"
        
    if not has_metal and not {'As', 'Sb'}.intersection(metal_ions):
        return False, "No essential mineral-forming elements found"
    
    # Check for water molecules
    water_pattern = Chem.MolFromSmiles("O")
    water_count = len(mol.GetSubstructMatches(water_pattern))
    
    # Check for common mineral patterns
    if c_count == 0:  # Simple inorganic minerals
        if len(metal_ions) > 0 and len(nonmetal_ions) > 0:
            return True, f"Simple inorganic mineral containing {', '.join(metal_ions)}"
            
    # Carbonate, sulfate, phosphate minerals
    if c_count <= 1 and {'O'}.issubset(nonmetal_ions):
        if len(charged_atoms) > 0 and len(metal_ions) > 0:
            return True, f"Ionic mineral containing {', '.join(metal_ions)}"
            
    # Hydrated minerals
    if water_count > 0:
        if c_count <= 1 and len(metal_ions) > 0:
            return True, f"Hydrated mineral containing {', '.join(metal_ions)}"
            
    # Special case for sulfides and similar minerals
    if {'S'}.issubset(nonmetal_ions) and len(metal_ions) > 0:
        if c_count == 0:
            return True, f"Sulfide mineral containing {', '.join(metal_ions)}"
            
    # Handle silicates
    if {'Si', 'O'}.issubset(nonmetal_ions) and len(metal_ions) > 0:
        if c_count == 0:
            return True, f"Silicate mineral containing {', '.join(metal_ions)}"
    
    return False, "Structure not consistent with mineral characteristics"