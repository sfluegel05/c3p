"""
Classifies: CHEBI:32952 amine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_amine(smiles: str):
    """
    Determines if a molecule is an amine based on its SMILES string.
    An amine is a compound formally derived from ammonia by replacing one, two, or three hydrogen atoms with hydrocarbyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Preprocess exclusion patterns
    nitro_pattern = Chem.MolFromSmarts('[N+](=O)[O-]')  # Nitro groups
    amide_pattern = Chem.MolFromSmarts('[C,c](=O)[N,n]')  # Amides
    sulfonamide_pattern = Chem.MolFromSmarts('[S,s](=O)(=O)[N,n]')  # Sulfonamides
    nitrile_pattern = Chem.MolFromSmarts('[C,c]#[N,n]')  # Nitriles
    azo_pattern = Chem.MolFromSmarts('[N,n]=[N,n]')  # Azo groups
    
    exclusion_matches = set()
    for pattern in [nitro_pattern, amide_pattern, sulfonamide_pattern, nitrile_pattern, azo_pattern]:
        if pattern:
            matches = mol.GetSubstructMatches(pattern)
            for match in matches:
                exclusion_matches.update(match)
    
    # Get all nitrogen atoms
    nitrogens = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]
    if not nitrogens:
        return False, "No nitrogen atoms present"
    
    # Add hydrogens for accurate valence calculation
    mol = Chem.AddHs(mol)
    
    for atom in nitrogens:
        idx = atom.GetIdx()
        if idx in exclusion_matches:
            continue  # Skip excluded functional groups
        
        # Check if part of aromatic system (e.g., pyridine, pyrrole)
        if atom.GetIsAromatic():
            continue
        
        # Check valence and charge conditions
        total_degree = atom.GetTotalDegree()
        formal_charge = atom.GetFormalCharge()
        
        # Valid amine conditions: either neutral with 3 bonds or +1 charge with 4 bonds
        if not ((formal_charge == 0 and total_degree == 3) or 
                (formal_charge == 1 and total_degree == 4)):
            continue
        
        # Count carbon and hydrogen substituents
        carbon_count = 0
        hydrogen_count = 0
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:
                carbon_count += 1
            elif neighbor.GetAtomicNum() == 1:
                hydrogen_count += 1
            else:  # Non-carbon/hydrogen substituent
                break
        else:  # Only executed if all neighbors are C or H
            if carbon_count < 1:
                continue  # At least one hydrocarbyl group required
            
            # Verify substituent count matches valence
            if formal_charge == 0 and (carbon_count + hydrogen_count) == 3:
                amine_type = ['primary', 'secondary', 'tertiary'][carbon_count - 1]
                return True, f"Contains {amine_type} amine group"
            elif formal_charge == 1 and (carbon_count + hydrogen_count) == 4:
                amine_type = ['primary', 'secondary', 'tertiary'][carbon_count - 1]
                return True, f"Contains {amine_type} ammonium group"
    
    return False, "No valid amine groups detected"