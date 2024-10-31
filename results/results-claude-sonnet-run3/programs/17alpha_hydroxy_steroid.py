from rdkit import Chem
from rdkit.Chem import AllChem

def is_17alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 17alpha-hydroxy steroid.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 17alpha-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic steroid core SMARTS pattern with specific stereochemistry at C17 and OH group
    steroid_pattern = '[#6]~1~2~[#6]~[#6]~[#6]~3~[#6](~[#6]~[#6]~[#6]~4~[#6]~[#6]~[#6]~[#6]~[#6]~4~[#6]~3~[#6]~[#6]~2)~[#6]~[#6]~[#6](@[H])(@[O])~[#6]~1'
    
    # Convert SMILES to 3D structure for better stereochemistry handling
    try:
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
    except:
        pass  # Continue even if 3D conversion fails

    # Check for the presence of the steroid core
    if not mol.HasSubstructMatch(Chem.MolFromSmarts(steroid_pattern)):
        return False, "Does not contain basic steroid core structure"

    # SMARTS pattern specifically for 17alpha-hydroxy group
    alpha_oh_pattern = '[C;R1]1([C;R1][C;R1][C;R1]2)([H])[C;R1]([H])[C;R1]([O;H1])[C;R1]2'
    
    if mol.HasSubstructMatch(Chem.MolFromSmarts(alpha_oh_pattern)):
        # Additional verification of stereochemistry
        matches = mol.GetSubstructMatches(Chem.MolFromSmarts(alpha_oh_pattern))
        for match in matches:
            c17_idx = match[2]  # Index of C17 in the pattern
            oh_idx = [i for i in match if mol.GetAtomWithIdx(i).GetSymbol() == 'O'][0]
            
            # Check if the OH is in alpha position
            bond = mol.GetBondBetweenAtoms(c17_idx, oh_idx)
            if bond is not None:
                return True, "Contains 17alpha-hydroxy steroid core with correct stereochemistry"
    
    return False, "Does not contain 17alpha-hydroxy steroid core structure"
# Pr=None
# Recall=0.0