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
    
    # Define common elements in minerals
    mineral_metals = {
        11: 'Na', 12: 'Mg', 13: 'Al', 19: 'K', 20: 'Ca', 
        26: 'Fe', 28: 'Ni', 29: 'Cu', 30: 'Zn', 56: 'Ba',
        55: 'Cs', 51: 'Sb', 33: 'As', 57: 'La', 46: 'Pd'
    }
    
    mineral_nonmetals = {
        8: 'O', 9: 'F', 15: 'P', 16: 'S', 17: 'Cl',
        14: 'Si', 5: 'B', 7: 'N'
    }
    
    # Track atom types and properties
    metal_atoms = []
    nonmetal_ions = set()
    c_atoms = []
    total_charge = 0
    water_count = 0
    
    # Analyze atoms
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        formal_charge = atom.GetFormalCharge()
        total_charge += formal_charge
        
        if atomic_num in mineral_metals:
            metal_atoms.append((atomic_num, formal_charge))
        elif atomic_num in mineral_nonmetals:
            nonmetal_ions.add(mineral_nonmetals[atomic_num])
        elif atomic_num == 6:  # Carbon
            c_atoms.append(atom)
            
    if len(metal_atoms) == 0 and not {'As', 'Sb'}.intersection(nonmetal_ions):
        return False, "No mineral-forming metal elements found"

    # Count waters of hydration
    water_pattern = Chem.MolFromSmiles("O")
    water_count = len(mol.GetSubstructMatches(water_pattern))
    
    # Analyze carbon environment
    organic_ligands = []
    inorganic_carbons = 0
    
    for c_atom in c_atoms:
        # Check if carbon is part of carbonate
        if (c_atom.GetDegree() == 3 and 
            sum(1 for n in c_atom.GetNeighbors() if n.GetAtomicNum() == 8) == 3):
            inorganic_carbons += 1
            continue
            
        # Check if carbon is part of legitimate organic ligand
        neighbors = c_atom.GetNeighbors()
        if any(n.GetAtomicNum() == 8 and 
              any(nn.GetAtomicNum() in mineral_metals.keys() for nn in n.GetNeighbors())
              for n in neighbors):
            organic_ligands.append(c_atom)
            continue
            
        return False, "Contains non-mineral organic components"

    # Validate charge balance
    if abs(total_charge) > 0.1:  # Allow for small floating point errors
        return False, "Unbalanced charges in structure"

    # Classify based on structural patterns
    if len(c_atoms) == 0:  # Simple inorganic minerals
        metal_names = [mineral_metals[num] for num, _ in metal_atoms]
        return True, f"Inorganic mineral containing {', '.join(set(metal_names))}"
        
    if inorganic_carbons > 0:  # Carbonates
        metal_names = [mineral_metals[num] for num, _ in metal_atoms]
        return True, f"Carbonate mineral containing {', '.join(set(metal_names))}"
        
    if len(organic_ligands) > 0:  # Minerals with organic ligands
        if all(c.GetDegree() <= 4 for c in organic_ligands):
            metal_names = [mineral_metals[num] for num, _ in metal_atoms]
            return True, f"Metal {metal_names[0]} mineral with organic ligands"
            
    if water_count > 0:  # Hydrated minerals
        if len(c_atoms) == 0 or inorganic_carbons > 0:
            metal_names = [mineral_metals[num] for num, _ in metal_atoms]
            return True, f"Hydrated mineral containing {', '.join(set(metal_names))}"
    
    return False, "Structure not consistent with mineral characteristics"