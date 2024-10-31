from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.AllChem import GetMorganFingerprintAsBitVect
from rdkit.DataStructs.cDataStructs import TanimotoSimilarity

def is_erythromycin(smiles: str):
    """
    Determines if a molecule is an erythromycin derivative.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an erythromycin, False otherwise
        str: Reason for classification
    """
    # Reference structure - erythromycin A
    erythromycin_a = "CC[C@H]1OC(=O)[C@H](C)[C@@H](O[C@H]2C[C@@](C)(OC)[C@@H](O)[C@H](C)O2)[C@H](C)[C@@H](O[C@@H]2O[C@H](C)C[C@@H]([C@H]2O)N(C)C)[C@](C)(O)C[C@@H](C)C(=O)[C@H](C)[C@@H](O)[C@]1(C)O"
    
    # Parse input molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"
        
    # Parse reference molecule
    ref_mol = Chem.MolFromSmiles(erythromycin_a)
    
    # Calculate Morgan fingerprints (ECFP4)
    query_fp = GetMorganFingerprintAsBitVect(mol, 2, 2048)
    ref_fp = GetMorganFingerprintAsBitVect(ref_mol, 2, 2048)
    
    # Calculate Tanimoto similarity
    similarity = TanimotoSimilarity(query_fp, ref_fp)
    
    # Check molecular formula
    query_formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
    ref_formula = Chem.rdMolDescriptors.CalcMolFormula(ref_mol)
    
    # Core characteristics of erythromycins:
    # - Large macrolide structure (similar to erythromycin A)
    # - Similar molecular formula to erythromycin A (C37H67NO13)
    # - High structural similarity to erythromycin A
    
    if similarity >= 0.85:
        return True, f"High structural similarity to erythromycin A (Tanimoto={similarity:.2f})"
    elif similarity >= 0.75:
        return True, f"Moderate structural similarity to erythromycin A (Tanimoto={similarity:.2f}), likely an erythromycin derivative"
    else:
        return False, f"Low structural similarity to erythromycin A (Tanimoto={similarity:.2f})"
# Pr=1.0
# Recall=0.6666666666666666