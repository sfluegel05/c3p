from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_bile_acid(smiles: str):
    """
    Determines if a molecule is a bile acid based on its structure.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a bile acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None
        
    # Check for steroid core (cyclopentanophenanthrene)
    steroid_core = False
    rings = mol.GetRingInfo()
    if rings.NumRings() >= 4:
        for ring in rings.AtomRings():
            if len(ring) == 6: # Check for 6-membered rings
                steroid_core = True
                
    if not steroid_core:
        return False, "Missing steroid core structure"
        
    # Check for carboxylic acid group (-COOH)
    carboxyl = False
    patt = Chem.MolFromSmarts('[CX3](=O)[OX2H1]')
    if mol.HasSubstructMatch(patt):
        carboxyl = True
        
    if not carboxyl:
        return False, "Missing carboxylic acid group"
        
    # Check for hydroxyl groups (-OH)
    hydroxyls = 0
    patt = Chem.MolFromSmarts('[OX2H1]')
    hydroxyls = len(mol.GetSubstructMatches(patt))
    
    if hydroxyls < 1:
        return False, "Missing hydroxyl groups"
        
    # Check for 5beta configuration
    # This is a simplification - full stereochemistry check would be more complex
    chiral_centers = rdMolDescriptors.CalcNumAtomStereoCenters(mol)
    if chiral_centers < 5:
        return False, "Incorrect stereochemistry"
        
    # Check molecular weight range typical for bile acids
    mw = Descriptors.ExactMolWt(mol)
    if not (350 < mw < 600):
        return False, "Molecular weight outside typical range for bile acids"
        
    return True, "Matches bile acid structural requirements: steroid core, carboxylic acid, hydroxyls, stereochemistry"
# Pr=0.9863481228668942
# Recall=0.6980676328502415