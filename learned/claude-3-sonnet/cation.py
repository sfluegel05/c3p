"""
Classifies: CHEBI:36916 cation
"""
"""
Classifies: CHEBI:36916 cation
"""
from rdkit import Chem
from rdkit.Chem import rdmolops

def is_cation(smiles: str):
    """
    Determines if a molecule is a cation based on its SMILES string.
    A cation is a monoatomic or polyatomic species having one or more 
    elementary charges of the proton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cation, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Calculate formal charges on each atom
    positive_charges = 0
    positive_atoms = []
    
    for atom in mol.GetAtoms():
        formal_charge = atom.GetFormalCharge()
        if formal_charge > 0:
            positive_charges += formal_charge
            positive_atoms.append(f"{atom.GetSymbol()}{atom.GetIdx()+1}")
    
    # Look for specific cationic patterns
    cationic_patterns = [
        Chem.MolFromSmarts('[N+]'), # Any positively charged nitrogen
        Chem.MolFromSmarts('[n+]'), # Any positively charged aromatic nitrogen
        Chem.MolFromSmarts('[P+]'), # Phosphonium
        Chem.MolFromSmarts('[S+]'), # Sulfonium
        Chem.MolFromSmarts('[O+]'), # Oxonium
        Chem.MolFromSmarts('[#7+]'), # Any positively charged nitrogen (alternative notation)
        Chem.MolFromSmarts('[NH4+]'), # Ammonium
        Chem.MolFromSmarts('[NH3+]'), # Primary ammonium
        Chem.MolFromSmarts('[NH2+]'), # Secondary ammonium
        Chem.MolFromSmarts('[NH+]'),  # Tertiary ammonium
        Chem.MolFromSmarts('[N+](C)(C)(C)'), # Quaternary ammonium
        Chem.MolFromSmarts('[Li+,Na+,K+,Rb+,Cs+,Fr+]'), # Alkali metals
        Chem.MolFromSmarts('[Be+2,Mg+2,Ca+2,Sr+2,Ba+2,Ra+2]'), # Alkaline earth metals
        Chem.MolFromSmarts('[#+1,#+2,#+3,#+4]') # Any atom with positive charge
    ]
    
    # If no explicit positive charges found, check for matches to cationic patterns
    if positive_charges == 0:
        for pattern in cationic_patterns:
            if pattern is not None and mol.HasSubstructMatch(pattern):
                matches = mol.GetSubstructMatches(pattern)
                for match in matches:
                    atom_idx = match[0]
                    atom = mol.GetAtomWithIdx(atom_idx)
                    positive_atoms.append(f"{atom.GetSymbol()}{atom_idx+1}")
                positive_charges += len(matches)
                
    # If still no positive charges found
    if positive_charges == 0:
        return False, "No positive charges or cationic groups found"
    
    # Success case - report all positive charges found
    charge_locations = ", ".join(positive_atoms)
    charge_str = f"{positive_charges}+" if positive_charges > 1 else "1+"
    return True, f"Found {charge_str} charge with cationic centers on: {charge_locations}"