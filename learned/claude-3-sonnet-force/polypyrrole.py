"""
Classifies: CHEBI:38077 polypyrrole
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_polypyrrole(smiles: str):
    """
    Determines if a molecule is a polypyrrole based on its SMILES string.
    A polypyrrole is a compound composed of two or more pyrrole units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polypyrrole, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 100:
        return False, "Molecular weight too low for polypyrrole"
    
    # Look for pyrrole rings
    pyrrole_pattern = Chem.MolFromSmarts("[nH]1cccc1")
    pyrrole_matches = mol.GetSubstructMatches(pyrrole_pattern)
    
    # Check for at least two pyrrole units
    if len(pyrrole_matches) < 2:
        return False, "Less than two pyrrole units found"
    
    # Check for conjugated system
    conjugated = False
    kekulization = AllChem.Kekulize(mol)
    for ring in kekulization.GetRingInfo().AtomRings():
        ring_atoms = [mol.GetAtomWithIdx(idx).GetSymbol() for idx in ring]
        if all(atom == 'C' or atom == 'N' for atom in ring_atoms):
            conjugated = True
            break
    
    if not conjugated:
        return False, "No conjugated ring system found"
    
    # Handle specific cases
    if "Co" in [atom.GetSymbol() for atom in mol.GetAtoms()]:
        return True, "Contains cobalt atom and pyrrole units (e.g., cobalt-precorrin-5B)"
    
    if any(smi in smiles for smi in ["C1=CC(=Cc2ccc[nH]2)N=C1", "c1ccc(Cc2ccc[nH]2)[nH]1"]):
        return True, "Contains dipyrrin or similar substructure"
    
    return True, "Contains two or more pyrrole units in a conjugated ring system"