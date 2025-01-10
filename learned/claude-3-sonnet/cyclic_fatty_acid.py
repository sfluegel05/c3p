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
    if rdMolDescriptors.CalcNumRings(mol) == 0:
        return False, "No rings found in structure"
        
    # Count atoms
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    num_nitrogens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    num_oxygens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    # Check for peptides (too many nitrogens)
    if num_nitrogens > 2:
        return False, "Too many nitrogens for a fatty acid"
    
    # Check for carbohydrates (too many oxygens relative to carbons)
    if num_oxygens > num_carbons/2:
        return False, "Too many oxygens relative to carbons"
    
    # Check for fatty acid characteristic groups
    carboxylic_acid = Chem.MolFromSmarts('[CX3](=O)[OX2H1]')
    ester = Chem.MolFromSmarts('[CX3](=O)[OX2][CX4]')
    
    has_acid = mol.HasSubstructMatch(carboxylic_acid)
    has_ester = mol.HasSubstructMatch(ester)
    
    if not (has_acid or has_ester):
        return False, "No carboxylic acid or ester group found"
    
    # Check for long carbon chain
    alkyl_chain = Chem.MolFromSmarts('CCCC')  # At least 4 consecutive carbons
    if not mol.HasSubstructMatch(alkyl_chain):
        return False, "No significant alkyl chain found"
    
    # Get ring information
    ring_info = mol.GetRingInfo()
    ring_sizes = [len(ring) for ring in ring_info.AtomRings()]
    
    # Fatty acids typically don't have very large rings
    if any(size > 8 for size in ring_sizes):
        return False, "Ring size too large for typical fatty acid"
    
    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 130 or mol_wt > 600:  # Upper limit to exclude complex compounds
        return False, "Molecular weight outside typical fatty acid range"
    
    # Compile ring description
    ring_desc = []
    for size in sorted(set(ring_sizes)):
        count = ring_sizes.count(size)
        ring_desc.append(f"{count} {size}-membered" + (" ring" if count == 1 else " rings"))
    
    ring_details = ", ".join(ring_desc)
    
    # Check ring types
    num_aromatic_rings = sum(1 for ring in ring_info.AtomRings() 
                            if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring))
    
    ring_type = "aromatic" if num_aromatic_rings > 0 else "aliphatic"
    if num_aromatic_rings > 0 and len(ring_sizes) > num_aromatic_rings:
        ring_type = "mixed aromatic and aliphatic"
    
    acid_type = "carboxylic acid" if has_acid else "ester"
    
    return True, f"Contains {ring_details} ({ring_type}) and a {acid_type} group"