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
    
    # Check for carboxylic acid group
    carboxylic_acid = Chem.MolFromSmarts('[CX3](=O)[OX2H1]')
    if not mol.HasSubstructMatch(carboxylic_acid):
        return False, "No carboxylic acid group found"
    
    # Check if molecule contains any rings
    if rdMolDescriptors.CalcNumRings(mol) == 0:
        return False, "No rings found in structure"
        
    # Get basic molecular properties
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Fatty acids typically have at least 8 carbons
    if num_carbons < 8:
        return False, "Too few carbons to be a fatty acid"
        
    # Check for reasonable molecular weight (typical fatty acids >130 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 130:
        return False, "Molecular weight too low for fatty acid"
    
    # Get ring information
    ring_info = mol.GetRingInfo()
    ring_sizes = [len(ring) for ring in ring_info.AtomRings()]
    
    # Compile ring information for the reason string
    ring_desc = []
    for size in sorted(set(ring_sizes)):
        count = ring_sizes.count(size)
        ring_desc.append(f"{count} {size}-membered" + (" ring" if count == 1 else " rings"))
    
    ring_details = ", ".join(ring_desc)
    
    # Additional check for aromatic rings
    num_aromatic_rings = sum(1 for ring in ring_info.AtomRings() 
                            if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring))
    
    ring_type = "aromatic" if num_aromatic_rings > 0 else "aliphatic"
    if num_aromatic_rings > 0 and len(ring_sizes) > num_aromatic_rings:
        ring_type = "mixed aromatic and aliphatic"
    
    return True, f"Contains {ring_details} ({ring_type}) and a carboxylic acid group"