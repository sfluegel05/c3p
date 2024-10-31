from rdkit import Chem
from rdkit.Chem import AllChem

def is_metal_carbonyl(smiles: str):
    """
    Determines if a molecule is a metal carbonyl complex.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a metal carbonyl, False otherwise
        str: Reason for classification
    """
    # Use non-sanitized molecule first to preserve metal bonds
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Try to sanitize, but keep going even if it fails
    try:
        Chem.SanitizeMol(mol, sanitizeOps=Chem.SANITIZE_ALL^Chem.SANITIZE_KEKULIZE)
    except:
        pass
        
    # List of metal atoms that commonly form carbonyls
    metals = {'Fe', 'Cr', 'Mo', 'W', 'Mn', 'Tc', 'Re', 'Ru', 'Os', 'Co', 'Rh', 'Ir', 'Ni', 'Pd', 'Pt', 'V'}
    
    # Find metal atoms in molecule
    metal_atoms = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() in metals:
            metal_atoms.append(atom)
            
    if not metal_atoms:
        return False, "No metal atoms found"
        
    # Check for CO ligands
    co_count = 0
    for metal in metal_atoms:
        for neighbor in metal.GetNeighbors():
            # Check if neighbor is carbon
            if neighbor.GetSymbol() == 'C':
                # Check if carbon has oxygen neighbor
                for c_neighbor in neighbor.GetNeighbors():
                    if c_neighbor.GetSymbol() == 'O':
                        # Get the bond between C and O
                        bond = mol.GetBondBetweenAtoms(neighbor.GetIdx(), c_neighbor.GetIdx())
                        # Check for various representations of CO bonds
                        if (bond.GetBondType() in [Chem.BondType.TRIPLE, Chem.BondType.DOUBLE] or 
                            (bond.GetBondType() == Chem.BondType.SINGLE and c_neighbor.GetFormalCharge() == 1)):
                            co_count += 1
                            break
                            
    if co_count == 0:
        return False, "No CO ligands found"
        
    metal_symbols = [atom.GetSymbol() for atom in metal_atoms]
    return True, f"Metal carbonyl complex with {len(metal_atoms)} {'/'.join(set(metal_symbols))} center(s) and {co_count} CO ligand(s)"
# Pr=1.0
# Recall=0.9444444444444444