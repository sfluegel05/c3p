from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize

def is_hydrazone(smiles: str):
    """
    Determines if a molecule contains a hydrazone group (R2C=NNR2).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule contains hydrazone group, False otherwise
        str: Reason for classification
    """
    # Check for valid SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"
        
    # Standardize the molecule (normalize tautomers etc)
    clean_mol = rdMolStandardize.Cleanup(mol)
    
    # SMARTS patterns for hydrazone group
    # Matches C=N-N where the carbons can be aliphatic or aromatic
    # and the terminal N can have various substituents
    hydrazone_patterns = [
        '[C,c]=[N]-[N]',  # General hydrazone pattern
        '[C,c]=[N]-[N]C(=[O,S])',  # Hydrazones with carbonyl/thiocarbonyl
        '[C,c]=[N]-[N]C#N',  # Cyano-hydrazones
        '[C,c]=[N]-[N]c1[c,n][c,n][c,n][c,n][c,n]1'  # Aromatic substituted hydrazones
    ]
    
    for pattern in hydrazone_patterns:
        substructure = Chem.MolFromSmarts(pattern)
        if substructure is None:
            continue
            
        matches = clean_mol.GetSubstructMatches(substructure)
        if matches:
            # Get atoms involved in match
            match_atoms = matches[0]
            # Get carbon atom that forms C=N bond
            carbon = clean_mol.GetAtomWithIdx(match_atoms[0])
            # Get first nitrogen (hydrazone N)
            n1 = clean_mol.GetAtomWithIdx(match_atoms[1])
            
            # Verify C=N double bond exists
            for bond in carbon.GetBonds():
                if (bond.GetOtherAtomIdx(carbon.GetIdx()) == n1.GetIdx() and 
                    bond.GetBondType() == Chem.BondType.DOUBLE):
                    substituents = []
                    # Check substituents on carbon
                    for neighbor in carbon.GetNeighbors():
                        if neighbor.GetIdx() != n1.GetIdx():
                            substituents.append(neighbor.GetSymbol())
                    
                    return True, f"Contains hydrazone group with C substituents: {','.join(substituents)}"
                    
    return False, "Does not contain hydrazone group"
# Pr=None
# Recall=None