from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdmolops import GetAdjacencyMatrix
import numpy as np

def is_labdane_diterpenoid(smiles: str):
    """
    Determines if a molecule is a labdane diterpenoid based on structural features.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a labdane diterpenoid, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Check for basic molecular properties
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')
    if num_carbons < 18 or num_carbons > 30:  # Allow for some flexibility due to substituents
        return False, "Number of carbons not consistent with labdane skeleton"

    # Generate 3D conformation for structure analysis
    try:
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
    except:
        pass  # Continue even if 3D conformation generation fails

    # Check for decalin core (bicyclic system)
    rings = mol.GetRingInfo()
    ring_atoms = rings.AtomRings()
    if len(ring_atoms) < 2:
        return False, "Missing required bicyclic ring system"

    # Find 6-membered rings
    six_membered_rings = [ring for ring in ring_atoms if len(ring) == 6]
    if len(six_membered_rings) < 2:
        return False, "Missing required two 6-membered rings"

    # Check for connected 6-membered rings (decalin system)
    if len(six_membered_rings) >= 2:
        adj_matrix = GetAdjacencyMatrix(mol)
        found_decalin = False
        for i, ring1 in enumerate(six_membered_rings):
            for ring2 in six_membered_rings[i+1:]:
                shared_atoms = set(ring1).intersection(set(ring2))
                if len(shared_atoms) == 2:  # Two rings share exactly 2 atoms
                    found_decalin = True
                    break
            if found_decalin:
                break
        if not found_decalin:
            return False, "Six-membered rings not properly fused"

    # Check for isoprenoid-like branching pattern
    features = []
    has_isoprenoid = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and len(atom.GetNeighbors()) == 4:
            methyl_count = sum(1 for n in atom.GetNeighbors() 
                             if n.GetSymbol() == 'C' and len(n.GetNeighbors()) == 1)
            if methyl_count >= 1:
                has_isoprenoid = True
                features.append("isoprenoid branching")
                break

    # Check for common functionalities
    if any(atom.GetSymbol() == 'O' for atom in mol.GetAtoms()):
        features.append("oxygen functionality")
    
    # Check for double bonds
    double_bonds = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            double_bonds += 1
    if double_bonds > 0:
        features.append(f"{double_bonds} double bond(s)")

    if features:
        feature_str = ", ".join(features)
        return True, f"Labdane diterpenoid with {feature_str}"
    else:
        return True, "Basic labdane diterpenoid skeleton"
# Pr=0.9285714285714286
# Recall=1.0