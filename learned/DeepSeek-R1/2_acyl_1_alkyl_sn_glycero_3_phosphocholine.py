"""
Classifies: CHEBI:36702 2-acyl-1-alkyl-sn-glycero-3-phosphocholine
"""
"""
Classifies: 2-acyl-1-alkyl-sn-glycero-3-phosphocholine (CHEBI:XXXXX)
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_2_acyl_1_alkyl_sn_glycero_3_phosphocholine(smiles: str):
    """
    Determines if a molecule is a 2-acyl-1-alkyl-sn-glycero-3-phosphocholine based on its SMILES string.
    Must have:
    - sn-glycerol backbone with stereochemistry
    - Phosphocholine group at the sn-3 position
    - Ether-linked alkyl chain at sn-1
    - Ester-linked acyl chain at sn-2
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Standardize to avoid stereochemistry issues
    Chem.SanitizeMol(mol)
    
    # Define core structure pattern using SMARTS
    # [1]: Phosphocholine group
    # [2]: Central chiral carbon (sn-2)
    # [3]: sn-1 ether oxygen
    # [4]: sn-2 ester oxygen
    core_pattern = Chem.MolFromSmarts("""
        [O-]P(=O)(OCC[N+](C)(C)C)O[C@H](CO[!H0])OC(=O)[!H0]
    """)
    
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Core structure mismatch"
    
    # Verify alkyl chain (ether) at sn-1
    ether_pattern = Chem.MolFromSmarts("[CX4][OX2][CX4][CH2X4]")  # Ether linkage
    if not mol.HasSubstructMatch(ether_pattern):
        return False, "Missing alkyl ether chain at sn-1"

    # Verify acyl chain (ester) at sn-2
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2][C@H]")  # Ester linkage at chiral center
    if not mol.HasSubstructMatches(ester_pattern):
        return False, "Missing acyl ester chain at sn-2"

    # Check molecular weight consistency (phosphocholine base ~250Da + chains)
    mol_wt = AllChem.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, "Molecular weight too low for typical structure"

    return True, "Has 1-alkyl-2-acyl-sn-glycero-3-phosphocholine structure"