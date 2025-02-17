"""
Classifies: CHEBI:17408 monoacylglycerol
"""
"""
Classifies: Monoacylglycerol
Definition: A glyceride in which any one of the R groups (position not specified) is an acyl group
while the remaining two R groups can be either H or alkyl groups.
In other words, a monoacylglycerol has a glycerol (three‐carbon) backbone where exactly one of
the hydroxyl groups is esterified while the other two remain as free hydroxyl groups (or are alkylated).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monoacylglycerol(smiles: str):
    """
    Determines if a molecule is a monoacylglycerol based on its SMILES string.
    
    Approach:
    1. Parse SMILES.
    2. Count the number of ester bonds by looking for the substructure –O–C(=O)–.
    3. Count free hydroxyl groups (–OH) present.
    4. Look for a characteristic glycerol backbone fragment.
       (Because one OH is replaced by the acyl group, a MAG should have evidence of a three‐carbon chain 
       with two free hydroxyl groups.)
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a monoacylglycerol, False otherwise.
        str: Explanation for the classification decision.
    """
    
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a pattern for an ester group: the oxygen connected to the acyl carbon.
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=O)")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    num_esters = len(ester_matches)
    
    if num_esters != 1:
        # If there are zero or more than one ester groups, it is not a monoacylglycerol.
        return False, f"Number of ester groups is {num_esters}; expected exactly 1 for a monoacylglycerol"
    
    # Count free hydroxyl groups: an -OH pattern.
    free_oh_pattern = Chem.MolFromSmarts("[OX2H]")
    free_oh_matches = mol.GetSubstructMatches(free_oh_pattern)
    num_free_oh = len(free_oh_matches)
    
    if num_free_oh < 2:
        # A glycerol derivative should have two free –OH groups
        return False, f"Found only {num_free_oh} free hydroxyl groups; expected at least 2 on a glycerol backbone"
    
    # Look for a glycerol backbone motif.
    # For glycerol one would expect a 3‐carbon chain that has oxygen substituents.
    # In unmodified glycerol, the scaffold is HO–CH2–CHOH–CH2–OH.
    # In a MAG one OH is acylated (i.e. replaced by –OC(=O)R) so we try two patterns:
    # One pattern for free glycerol and one pattern for MAG backbone with one acyl group.
    glycerol_free_smarts = "OCC(O)CO"  # matches glycerol if all three OH were free
    glycerol_acyl_smarts = "OCC(OC(=O))CO"  # one OH replaced by ester (acyl) group
    glycerol_free = Chem.MolFromSmarts(glycerol_free_smarts)
    glycerol_acyl = Chem.MolFromSmarts(glycerol_acyl_smarts)
    
    has_backbone = mol.HasSubstructMatch(glycerol_free) or mol.HasSubstructMatch(glycerol_acyl)
    if not has_backbone:
        # If we cannot identify a glycerol-like backbone pattern then we cannot confidently assign MAG.
        return False, "Glycerol backbone pattern not found"
    
    # Optionally, one could add additional checks such as limits on molecular weight or count of rotatable bonds.
    # For many MAGs the weight is greater than 200 Da.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is too low for a monoacylglycerol"
    
    # Passed all heuristic checks.
    return True, "Contains a glycerol backbone with exactly one acyl (ester) group and two free hydroxyls, consistent with a monoacylglycerol"

# For testing purposes, one might call the function with one of the sample MAG SMILES strings.
if __name__ == "__main__":
    test_smiles = "CCCCCCCC(=O)OC[C@@H](O)CO"  # 1-octanoyl-sn-glycerol
    result, reason = is_monoacylglycerol(test_smiles)
    print(f"Result: {result}\nReason: {reason}")