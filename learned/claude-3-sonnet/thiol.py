"""
Classifies: CHEBI:29256 thiol
"""
"""
Classifies: CHEBI:33839 thiol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_thiol(smiles: str):
    """
    Determines if a molecule is a thiol based on its SMILES string.
    A thiol is an organosulfur compound with a -SH group attached to an aliphatic or aromatic carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a thiol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Generate tautomers
    tautomers = AllChem.EnumerateTautomers(mol)
    
    # Check each tautomer for thiol functional group
    for tautomer in tautomers:
        # Look for sulfhydryl (-SH) group
        sh_pattern = Chem.MolFromSmarts("[SH]")
        sh_matches = tautomer.GetSubstructMatches(sh_pattern)
        if not sh_matches:
            continue
        
        # Check if sulfhydryl is attached to aliphatic or aromatic carbon
        is_thiol = False
        for atom_idx in sh_matches[0]:
            atom = tautomer.GetAtomWithIdx(atom_idx)
            if atom.GetTotalNumHs() == 1:  # Check for -SH
                neighbors = atom.GetNeighbors()
                if any(nb.GetAtomicNum() == 6 and nb.GetHybridization() in (Chem.HybridizationType.SP2, Chem.HybridizationType.SP3) for nb in neighbors):
                    is_thiol = True
                    break
        
        if is_thiol:
            # Check for additional descriptors to improve classification
            mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
            if mol_wt < 50:  # Thiols typically have higher molecular weight
                continue
            num_sulfur = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16)
            if num_sulfur > 1:  # Avoid molecules with disulfide bridges
                continue
            
            return True, "Contains sulfhydryl (-SH) group attached to aliphatic or aromatic carbon"
    
    return False, "No thiol functional group found"