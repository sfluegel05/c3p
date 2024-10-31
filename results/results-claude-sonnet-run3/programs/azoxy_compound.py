from rdkit import Chem
from rdkit.Chem import AllChem

def is_azoxy_compound(smiles: str):
    """
    Determines if a molecule is an azoxy compound (RN=N(+)(O(-))R).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an azoxy compound, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for N=N(+)(O(-)) pattern using SMARTS
    azoxy_pattern = Chem.MolFromSmarts('[N]=[N+]([O-])')
    if mol.HasSubstructMatch(azoxy_pattern):
        # Find all matching atoms
        matches = mol.GetSubstructMatches(azoxy_pattern)
        
        for match in matches:
            n1, n2, o = match
            
            # Get the atoms
            n1_atom = mol.GetAtomWithIdx(n1)
            n2_atom = mol.GetAtomWithIdx(n2) 
            o_atom = mol.GetAtomWithIdx(o)
            
            # Verify formal charges
            if n2_atom.GetFormalCharge() != 1 or o_atom.GetFormalCharge() != -1:
                continue
                
            # Verify N=N double bond
            bond = mol.GetBondBetweenAtoms(n1, n2)
            if bond.GetBondType() != Chem.BondType.DOUBLE:
                continue
                
            # Verify N-O single bond
            bond = mol.GetBondBetweenAtoms(n2, o)
            if bond.GetBondType() != Chem.BondType.SINGLE:
                continue
                
            # Verify nitrogen atoms
            if n1_atom.GetSymbol() != 'N' or n2_atom.GetSymbol() != 'N':
                continue
                
            # Verify oxygen atom
            if o_atom.GetSymbol() != 'O':
                continue
                
            return True, "Contains N=N(+)(O(-)) azoxy group"
            
    return False, "Does not contain azoxy group"
# Pr=None
# Recall=None