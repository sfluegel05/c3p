"""
Classifies: CHEBI:59238 cyclic fatty acid
"""
"""
Classifies: cyclic fatty acid
Definition: Any fatty acid containing anywhere in its structure a ring of atoms.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_cyclic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a cyclic fatty acid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a cyclic fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for rings
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    if num_rings == 0:
        return False, "No rings found in structure"
    if num_rings > 5:
        return False, "Too many rings for a typical fatty acid"
        
    # Count atoms
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    num_nitrogens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    num_oxygens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    # Check atom ratios
    if num_nitrogens > 1:
        return False, "Too many nitrogens for a fatty acid"
    if num_oxygens > 4:
        return False, "Too many oxygens for a typical fatty acid"
    if num_carbons < 8:
        return False, "Too few carbons for a fatty acid"
    
    # Check for carboxylic acid group
    carboxylic_acid = Chem.MolFromSmarts('[CX3](=O)[OX2H1]')
    if not mol.HasSubstructMatch(carboxylic_acid):
        return False, "No carboxylic acid group found"
    
    # Get ring information
    ring_info = mol.GetRingInfo()
    ring_sizes = [len(ring) for ring in ring_info.AtomRings()]
    
    # Check ring sizes
    if any(size > 8 for size in ring_sizes):
        return False, "Ring size too large for typical fatty acid"
    if any(size < 3 for size in ring_sizes):
        return False, "Ring size too small to be stable"
    
    # Count aromatic rings
    num_aromatic_rings = sum(1 for ring in ring_info.AtomRings() 
                            if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring))
    
    if num_aromatic_rings > 1:
        return False, "Too many aromatic rings for a typical fatty acid"
    
    # Check for common non-fatty acid structures
    steroid_pattern = Chem.MolFromSmarts('C1CC2CCC3C(C2)CCC4C3(CCC4)C1')
    if mol.HasSubstructMatch(steroid_pattern):
        return False, "Steroid-like structure detected"
        
    # Compile ring description
    ring_desc = []
    for size in sorted(set(ring_sizes)):
        count = ring_sizes.count(size)
        ring_desc.append(f"{count} {size}-membered" + (" ring" if count == 1 else " rings"))
    
    ring_details = ", ".join(ring_desc)
    
    # Determine ring type
    ring_type = "aromatic" if num_aromatic_rings == len(ring_sizes) else "aliphatic"
    if num_aromatic_rings > 0 and num_aromatic_rings < len(ring_sizes):
        ring_type = "mixed aromatic and aliphatic"
    
    return True, f"Contains {ring_details} ({ring_type}) and a carboxylic acid group"