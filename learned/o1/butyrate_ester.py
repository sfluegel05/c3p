"""
Classifies: CHEBI:50477 butyrate ester
"""
"""
Classifies: butyrate ester

Definition: Any carboxylic ester where the carboxylic acid component is butyric acid.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_butyrate_ester(smiles: str):
    """
    Determines if a molecule is a butyrate ester based on its SMILES string.
    A butyrate ester is any carboxylic ester where the carboxylic acid component is butyric acid.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a butyrate ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for ester functional group
    ester_pattern = Chem.MolFromSmarts("C(=O)O[C;!$([C](=O))]")
    if ester_pattern is None:
        return False, "Failed to create ester pattern"
    
    # Find all ester groups in the molecule
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester groups found"
    
    # Generate butyric acid molecule for comparison
    butyric_acid_smiles = "CCCC(=O)O"
    butyric_acid = Chem.MolFromSmiles(butyric_acid_smiles)
    butyric_acid_frag = Chem.RemoveHs(butyric_acid)
    
    # Loop through each ester group and analyze the acyl component
    for match in ester_matches:
        ester_atoms = [mol.GetAtomWithIdx(idx) for idx in match]
        
        # Break the ester bond to get acyl and alkyl fragments
        editable_mol = Chem.EditableMol(mol)
        bond_to_break = mol.GetBondBetweenAtoms(match[1], match[2])
        if bond_to_break is None:
            continue
        editable_mol.RemoveBond(match[1], match[2])
        fragments = Chem.GetMolFrags(editable_mol.GetMol(), asMols=True)
        
        # Identify the acyl fragment (should contain the carbonyl carbon)
        acyl_fragment = None
        for frag in fragments:
            if frag.HasSubstructMatch(Chem.MolFromSmarts("C=O")):
                acyl_fragment = frag
                break
        
        if acyl_fragment is None:
            continue
        
        # Compare acyl fragment to butyric acid (without hydroxyl)
        acyl_smiles = Chem.MolToSmiles(acyl_fragment, isomericSmiles=True)
        # Add hydroxyl to make it comparable to butyric acid
        acyl_smiles += "O"
        acyl_mol = Chem.MolFromSmiles(acyl_smiles)
        if acyl_mol is None:
            continue
        # Use RDKit's Isomeric mol matching for comparison
        if acyl_mol.HasSubstructMatch(butyric_acid_frag):
            return True, "Contains butyrate ester group"
    
    return False, "Does not contain butyrate ester group"