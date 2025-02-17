"""
Classifies: CHEBI:16337 phosphatidic acid
"""
"""
Classifies: Phosphatidic acid
Definition: A derivative of glycerol in which one hydroxy group is esterified with phosphoric acid
and the other two hydroxy groups are esterified with fatty acids.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidic_acid(smiles: str):
    """
    Determines if a molecule is a phosphatidic acid based on its SMILES string.
    A phosphatidic acid contains one phosphate group (P bound to three oxygens),
    two fatty acid ester groups (OC(=O) fragments) attached to a glycerol backbone.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is phosphatidic acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Criterion 1: Check for the phosphate group.
    # In phosphatidic acids the phosphorus should be bonded to at least three oxygens.
    phosphate_found = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "P":
            oxy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetSymbol() == "O"]
            if len(oxy_neighbors) >= 3:
                phosphate_found = True
                break
    if not phosphate_found:
        return False, "Phosphate group not found (no P atom with at least 3 O neighbors)"
    
    # Criterion 2: Count fatty acid ester groups.
    # These are represented by the substructure "OC(=O)". This pattern should appear exactly twice.
    fatty_acid_pattern = Chem.MolFromSmarts("OC(=O)")
    fa_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fa_matches) != 2:
        return False, f"Found {len(fa_matches)} fatty acid ester group(s); exactly 2 are expected"
    
    # Criterion 3: Check for a glycerol backbone linkage.
    # One heuristic is to look for a fragment where a phosphate is linked through oxygen to 
    # a small carbon chain bearing the ester functions. For example, many phosphatidic acids 
    # show a pattern similar to "COC(=O)" in the branch away from the phosphate.
    # Here we adopt a heuristic SMARTS that attempts to capture a common glycerol-phosphate fragment:
    glycerol_pattern = Chem.MolFromSmarts("OCCOP(=O)(O)O")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "Glycerol backbone with phosphate linkage not found"
    
    # Additional sanity check: Many phosphatidic acids have a molecular weight above 300 Da.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, f"Molecular weight too low for phosphatidic acid ({mol_wt:.1f} Da)"
    
    return True, "Contains glycerol backbone with one phosphate and exactly two fatty acid chains"

# Example usage (you can remove or modify these lines for your integration):
if __name__ == "__main__":
    test_smiles = "P(OC[C@H](OC(=O)CCCCCCCCCCCCC)COC(=O)CCCCCCCCC/C=C\\CCCCCCCC)(O)(O)=O"
    result, reason = is_phosphatidic_acid(test_smiles)
    print("Is phosphatidic acid?", result)
    print("Reason:", reason)